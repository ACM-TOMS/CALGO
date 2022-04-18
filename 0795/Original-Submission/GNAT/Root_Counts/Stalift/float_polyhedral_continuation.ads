with text_io;                            use text_io;
with Integer_Vectors,Float_Vectors;
with Complex_Vectors_of_Vectors;
with Arrays_of_Float_Vector_Lists;       use Arrays_of_Float_Vector_Lists;
with Float_Mixed_Subdivisions;           use Float_Mixed_Subdivisions;
with Complex_Laurent_Polynomial_Systems; use Complex_Laurent_Polynomial_Systems;
with Laurent_Jacobi_Matrices;            use Laurent_Jacobi_Matrices;
with Exponent_Vectors;                   use Exponent_Vectors;
with Solutions;                          use Solutions;

package Float_Polyhedral_Continuation is

-- DESCRIPTION :
--   This package implements polyhedral homotopy continuation methods,
--   based on mixed subdivision induced by floating-point lifting.
--   The continuation is organized in three layers:
--     1. inner normal, tracking of paths for one poly
--     2. mixed cell, recursion is needed when the cell is not fine;
--     3. mixed subdivision, for all cells in the subdivision.
--   Each layer has two versions: a silent and a reporting version.

-- FIRST LAYER : polyhedral continuation for one transformation.

  procedure Mixed_Continuation
                ( mix : in Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  normal : in Float_Vectors.Vector;
                  sols : in out Solution_List );

  procedure Mixed_Continuation
                ( file : in file_type; mix : in Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  normal : in Float_Vectors.Vector;
                  sols : in out Solution_List );

  -- DESCRIPTION : polyhedral continuation with coefficient homotopy.

  -- ON ENTRY :
  --   file       file to write intermediate results on;
  --   mix        type of mixture;
  --   lifted     lifted supports of polynomial system, in original order;
  --   h          coefficient homotopy;
  --   c          coefficients of homotopy;
  --   j          coefficient Jacobi matrix;
  --   m          multiplication factors in coefficient Jacobi matrix;
  --   normal     normal to a mixed cell;
  --   sols       start solutions of the subsystem which corresponds
  --              with the mixed cell with given inner normal.

  -- ON RETURN :
  --   sols       the solutions of p, which correspond to one mixed cell.

-- SECOND LAYER : polyhedral continuaton for one mixed cell.

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List );

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral coefficient-homotopy continuation for one mixed cell.

  -- REQUIRED : polynomials in p must be ordered according to mix.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          non-lifted Laurent polynomial system;
  --   lifted     lifted supports, with original order of points;
  --   h          coefficient homotopy;
  --   c          coefficients of homotopy;
  --   e          the exponent vectors of the unlifted system;
  --   j          coefficient Jacobi matrix;
  --   m          multiplication factors in coefficient Jacobi matrix;
  --   mix        type of mixture;
  --   mic        a mixed cell.

  -- ON RETURN :
  --   sols       the solution list of p;
  --   sols_last  pointer to last element of the list sols.

-- THIRD LAYER : polyhedral continuation for a mixed subdivision.

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List );

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List );

  -- DESCRIPTION :
  --   Polyhedral coefficient-homotopy continuation for a mixed subdivision.

  -- REQUIRED : polynomials in p must be ordered according to mix.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          non-lifted Laurent polynomial system;
  --   lifted     lifted supports, in original order;
  --   h          coefficient homotopy;
  --   c          coefficients of homotopy;
  --   e          the exponent vectors of the unlifted system;
  --   j          coefficient Jacobi matrix;
  --   m          multiplication factors in coefficient Jacobi matrix;
  --   mix        type of mixture;
  --   mixsub     a collection of mixed cells.

  -- ON RETURN :
  --   sols       the solution list of p.

end Float_Polyhedral_Continuation;
