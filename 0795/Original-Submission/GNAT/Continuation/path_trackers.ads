with text_io;                          use text_io;
with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Float_Vectors;
with Complex_Numbers,Complex_Vectors;  use Complex_Numbers,Complex_Vectors;
with Complex_Matrices,Solutions;       use Complex_Matrices,Solutions;
with Continuation_Data;                use Continuation_Data;

package Path_Trackers is

-- DESCRIPTION :
--   This package offers some routines for tracking solution paths,
--   by using an increment-and-fix predictor-corrector method.

--   The following options can be made :
--    (Linear,Circular)
--       A linear path tracker takes t from a starting to a target value.
--       For computing winding numbers, a circular path tracker is needed.
--    (Single,Multiple)
--       A single path tracker only deals with one path at a time.
--       A multiple path tracker follows more than one path when it is called.
--    (Normal,Conditioned)
--       A normal path tracker does not compute an estimate for the inverse of
--       the condition number of the Jacobian matrix.  This additional work
--       is done by a conditioned path tracker.
--    (Silent,Reporting)
--       A silent path tracker does not produce any output on file.
--       A reporting path tracker allows to put intermediate results on file.

--   By combining these options, the following path trackers are provided:

--     Linear_Single_Normal_Silent_Continue
--     Linear_Single_Normal_Reporting_Continue
--     Linear_Single_Conditioned_Silent_Continue
--     Linear_Single_Conditioned_Reporting_Continue
--     Linear_Multiple_Normal_Silent_Continue
--     Linear_Multiple_Normal_Reporting_Continue
--     Linear_Multiple_Conditioned_Silent_Continue
--     Linear_Multiple_Conditioned_Reporting_Continue
--     Circular_Single_Normal_Reporting_Continue
--     Circular_Single_Conditioned_Reporting_Continue

-- LINEAR PATH FOLLOWING FOR ONE PATH :

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Single_Normal_Silent_Continue
                ( s : in out Solu_Info; target : in double_complex;
                  tol : in double_float; proj : in boolean;
                  p : in Pred_Pars; c : in Corr_Pars );
 
  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Single_Normal_Reporting_Continue
                ( file : in file_type; s : in out Solu_Info;
                  target : in double_complex; tol : in double_float;
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Single_Conditioned_Silent_Continue
                ( s : in out Solu_Info; target : in double_complex;
                  tol : in double_float; proj : in boolean;
                  rtoric : in natural;
                  v : in out Float_Vectors.Link_to_Vector;
                  errorv : in out double_float;
                  p : in Pred_Pars; c : in Corr_Pars );

  generic

    with function Norm ( x : Vector) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Single_Conditioned_Reporting_Continue
                ( file : in file_type; s : in out Solu_Info;
                  target : in double_complex; tol : in double_float;
                  proj : in boolean; rtoric : in natural;
                  v : in out Float_Vectors.Link_to_Vector;
                  errorv : in out double_float;
                  p : in Pred_Pars; c : in Corr_Pars );

  -- DESCRIPTION :
  --   This routine follows a path of solutions of the system
  --   H(x,t) = 0, with t : t ---> target.
  --   An increment and fix path following technique is applied.

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   s          start solution and initial value of t;
  --   target     target value of the continuation parameter;
  --   tol        tolerance to decide when two double_floats are the same;
  --   proj       when perpendicular-projective corrector has to be used;
  --   rtoric     order of extrapolation for computation of path directions;
  --   v          direction of toric compactificiation, null when (rtoric = 0);
  --   errorv     error on the current direction;
  --   p          parameters for the predictor;
  --   c          parameters for the corrector.

  -- ON RETURN :
  --   s          the computed solution of H(x,t) = 0;
  --   v          direction of the compactification, when rtoric > 0;
  --   errorv     difference with previously computed direction.

-- LINEAR PATH FOLLOWING FOR A NUMBER OF PATHS :

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Multiple_Normal_Silent_Continue
                ( s : in out Solu_Info_Array;
                  target : in double_complex; tol,dist_sols : in double_float;
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Multiple_Normal_Reporting_Continue
                ( file : in file_type; s : in out Solu_Info_Array;
                  target : in double_complex; tol,dist_sols : in double_float;
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Multiple_Conditioned_Silent_Continue
                ( s : in out Solu_Info_Array;
                  target : in double_complex; tol,dist_sols : in double_float;
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Linear_Multiple_Conditioned_Reporting_Continue
                ( file : in file_type; s : in out Solu_Info_Array;
                  target : in double_complex; tol,dist_sols : in double_float;
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );

  -- DESCRIPTION :
  --   This routine follows simultaneously a number of paths in order
  --   to avoid clustering of solutions.

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   s          array of start solutions, all for the same t;
  --   target     target value of the continuation parameters;
  --   tol        tolerance to decide when two double_floats are the same;
  --   dist_sols  distance to be kept between the solutions;
  --   proj       indicates whether perpendicular-projective corrector has
  --              to be used or not;
  --   p          parameters for the predictor;
  --   c          parameters for the corrector.

  -- ON RETURN :
  --   s          the computed solutions of H(x,t) = 0.

-- CIRCULAR PATH FOLLOWING FOR ONE PATH :

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Circular_Single_Normal_Reporting_Continue
                ( file : in file_type; s : in out Solu_Info;
                  target : in double_complex; tol,epslop : in double_float;
                  wc : out natural; max_wc : in natural;
                  sum,all_sum : out Vector;
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function H  ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Vector;
    with function dH ( x : Vector; t : double_complex ) return Matrix;

  procedure Circular_Single_Conditioned_Reporting_Continue
                ( file : in file_type; s : in out Solu_Info;
                  target : in double_complex; tol,epslop : in double_float;
                  wc : out natural; max_wc : in natural;
                  sum,all_sum : out Vector; 
                  proj : in boolean; p : in Pred_Pars; c : in Corr_Pars );


  -- DESCRIPTION :
  --   This routine follows a path of solutions of the system
  --   H(x,t) = 0, with t circuling around the target.
  --   An increment and fix path following technique is applied.

  -- ON ENTRY :
  --   file       to write intermediate results on;
  --   s          array of start solutions, all for the same t;
  --   target     target value of the continuation parameters;
  --   tol        tolerance to decide when two floats are the same;
  --   epslop     tolerance to decide when w(0) = w(2*PI*wc);
  --   max_wc     maximum bound for winding number;
  --   proj       indicates whether perpendicular-projective corrector has
  --              to be used or not;
  --   p          parameters for the predictor;
  --   c          parameters for the corrector.

  -- ON RETURN :
  --   s          the computed solution of H(x,t) = 0;
  --   wc         estimated cycle number;
  --   sum        is the trapezium sum over the equidistant points;
  --   all_sum    is the trapezium sum over all the points.

end Path_Trackers;
