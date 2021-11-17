with text_io;                           use text_io;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Float_Vectors;
with Float_Vectors_of_Vectors;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Matrices,Solutions;        use Complex_Matrices,Solutions;

package Increment_and_Fix_Continuation is

-- DESCRIPTION :
--   This package provides a general implementation of an increment-and-fix
--   continuation method.  The generic parameters are a norm, an evaluator and
--   a differentiator of the homotopy.
--   There are two basic versions: a silent and a reporting one.
--   The silent continuation simply performs its calculations without output
--   of intermediate results.  The reporting continuation routine allows to
--   put various kinds of intermediate results on a file.
--   It is assumed that the continuation parameters are already determined
--   before calling these routines (see Continuation_Parameters).
--   For both the silent and the reporting version, the facility is added
--   to estimate the directions of the solution paths, useful in a polyhedral
--   end game.

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Silent_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 target : in double_complex := CMPLX(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Reporting_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean;
                 target : in double_complex := CMPLX(1.0) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   proj      for projective-perpendicular path following;
  --   target    value for the continuation parameter at the end.
 
  -- ON RETURN :
  --   sols      the computed solutions.

-- WITH THE ESTIMATION OF THE PATH DIRECTIONS :

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Silent_Toric_Continue
               ( sols : in out Solution_List; proj : in boolean;
                 v : in out Float_Vectors_of_Vectors.Vector;
                 errv : in out Float_Vectors.Vector;
                 target : in double_complex := CMPLX(1.0) );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Reporting_Toric_Continue
               ( file : in file_type; sols : in out Solution_List;
                 proj : in boolean;
                 v : in out Float_Vectors_of_Vectors.Vector;
                 errv : in out Float_Vectors.Vector;
                 target : in double_complex := CMPLX(1.0) );

  -- DESCRIPTION :
  --   This routine implements the continuation strategy with the estimation
  --   of the directions of the solution paths at the end.

  -- ON ENTRY :
  --   file      to write intermediate results on (if Reporting_);
  --   sols      the start solutions;
  --   proj      for projective-perpendicular path following;
  --   v         v must be initialized with zero vectors
  --             and v'range is 1..Length_Of(sols);
  --   errv      errors on the computed directions;
  --   target    value for the continuation parameter at the end.

  -- ON RETURN :
  --   sols      the computed solutions;
  --   v         directions of the solution paths;
  --   errv      errors on the computed directions.

end Increment_and_Fix_Continuation;
