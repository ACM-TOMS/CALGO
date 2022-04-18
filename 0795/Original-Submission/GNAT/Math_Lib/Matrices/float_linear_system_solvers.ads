with Floating_Point_Numbers;          use Floating_Point_Numbers;
with Integer_Vectors,Float_Vectors;
with Float_Matrices;                  use Float_Matrices;

package Float_Linear_System_Solvers is

-- DESCRIPTION :
--   This package offers a few routines to solve linear systems of equations.
--   The code for lufac, lufco and lusolve is a literal translation from the
--   f77-linpack code.  We distinguish between static triangulators, which
--   take the matrix as a whole at one, and dynamic triangulators, which
--   allow to triangulate row after row.

-- STATIC TRIANGULATORS :

  procedure lufac ( a : in out Matrix; n : in integer;
                    ipvt : out Integer_Vectors.Vector; info : out integer );

  -- DESCRIPTION :
  --   lufac factors a float matrix by gaussian elimination

  --   lufac is usually called by lufco, but it can be called
  --   directly with a saving of time if rcond is not needed.
  --   (time for lufco) = (1 + 9/n)*(time for lufac).

  -- ON ENTRY :
  --   a       a floating point matrix(1..n,1..n) to be factored
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
  --   lufco factors a floating point matrix by gaussian elimination
  --   and estimates the condition of the matrix.

  -- If rcond is not needed, lufac is slightly faster.
  -- To solve a*x = b, follow lufco by lusolve.

  -- ON ENTRY :
  --   a       a floating point matrix(1..n,1..n) to be factored
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
                      b : in out Float_Vectors.Vector );

  -- DESCRIPTION :
  --   lusolve solves the system a*x = b using the factors
  --   computed by lufac or lufco

  -- ON ENTRY :
  --   a       a floating-point matrix(1..n,1..n), the output from
  --           lufac or lufco;
  --   n       the dimension of the matrix a;
  --   ipvt    the pivot vector from lufac or lufco;
  --   b       the right hand side vector.

  -- ON RETURN :
  --   b       the solution vector x.

  procedure Triangulate ( a : in out Matrix; n,m : in integer );

  -- DESCRIPTION :
  --   triangulate makes the n*m matrix a triangular using
  --   Gaussian elimination.

  -- ON ENTRY :
  --   a       a floating-point matrix(1..n,1..m);
  --   n       the number of rows of a;
  --   m       the number of columns of a.

  -- ON RETURN :
  --   a       the triangulated matrix.

  procedure Diagonalize ( a : in out Matrix; n,m : in integer );

  -- DESCRIPTION :
  --   Diagonalize makes the n*m floating-point matrix a diagonal.

  -- ON ENTRY :
  --   a       a floating-point matrix(1..n,1..m);
  --   n       the number of rows of a;
  --   m       the number of columns of a.

  -- ON RETURN :
  --   a       the diagonalized matrix.

-- DYNAMIC TRIANGULATORS :

  procedure Upper_Triangulate
               ( row : in natural; mat : in out Matrix; tol : in double_float;
                 ipvt : in out Integer_Vectors.Vector; pivot : out integer );

  -- DESCRIPTION :
  --   Makes the matrix upper triangular by updating the current row of
  --   the matrix.  If pivot = 0 on return, then the matrix is singular.

  -- REQUIRED :
  --   The matrix is upper triangular up to current row, which means that
  --   abs(max(i,i)) > tol, for i in mat'first(1)..row-1.
  --   The matrix might have more columns than rows, but ipvt'range should
  --   match the range of the first columns in the matrix.

  procedure Upper_Triangulate
               ( roweli : in natural; elim : in Matrix; tol : in double_float;
                 rowmat : in natural; mat : in out Matrix );
  procedure Upper_Triangulate
               ( roweli : in natural; elim : in Matrix; tol : in double_float;
                 firstrow,lastrow : in natural; mat : in out Matrix );

  -- DESCRIPTION :
  --   Using the matrix elim, the unknown roweli is eliminated in mat.

  -- ON ENTRY :
  --   roweli    current unknown to be eliminated;
  --   elim      elim(1..roweli,m'range(2)) is upper triangular;
  --   firstrow  indicates start in range of mat to be updated;
  --   lastrow   indicates end in range of mat to be updated;
  --   roweli    indicates the current unknown to be updated;
  --   mat       the unknows before roweli are already eliminated.

  -- ON RETURN :
  --   mat       the updated matrix after elimination of roweli.

  procedure Switch ( ipvt : in Integer_Vectors.Vector; row : in integer;
                     mat : in out Matrix );

  -- DESCRIPTION :
  --   Applies the pivoting information ipvt to the row of the matrix.

  procedure Switch ( k,pivot,first,last : in integer; mat : in out Matrix );

  -- DESCRIPTION :
  --   Swaps the column k with the pivot column for row first..last in mat.

  function Solve ( mat : Matrix; tol : double_float;
                   ipvt : Integer_Vectors.Vector ) return Float_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns a solution to the homogenous linear system with inequalities
  --   in the rows of the matrix.

  -- REQUIRED : the matrix is upper triangular.

end Float_Linear_System_Solvers;
