with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers,Complex_Vectors;
with Complex_Matrices,Solutions;         use Complex_Matrices,Solutions;
with Continuation_Data;                  use Continuation_Data;

package Predictors is

-- DESCRIPTION :
--   This package contains several implementations for the predictor 
--   in an increment-and-fix continuation.

--   The predictor provides a prediction both for the continuation parameter t
--   and for the solution(s) x.

--   For the continuation paramter t the following options can be made :
--     Real      : linear prediction, simply adds the step size;
--     Complex   : can make predictions in complex space;
--     Circular  : to perform a circular sample, for winding numbers;
--     Geometric : distances to target form geometric series.

--   For the solution vector x the following options are provided :
--     Secant  : linear extrapolation using differences;
--     Tangent : linear extrapolation using the first derivatives;
--     Hermite : third-order extrapolation using first derivatives.
--   Furthermore, these predictors for x can be applied
--   for one solution (Single) or for an array of solutions (Multiple).

--   By combining these options, the following 13 predictors are provided :

--     Secant_Single_Real_Predictor
--     Secant_Single_Complex_Predictor
--     Secant_Multiple_Real_Predictor
--     Secant_Multiple_Complex_Predictor
--     Tangent_Single_Real_Predictor
--     Tangent_Single_Complex_Predictor
--     Tangent_Multiple_Real_Predictor
--     Tangent_Multiple_Complex_Predictor

--     Secant_Circular_Predictor
--     Secant_Geometric_Predictor
--     Tangent_Circular_Predictor
--     Tangent_Geometric_Predictor

--     Hermite_Single_Real_Predictor

-- The order in which these predictors are listed depends on their mutual
-- resemblances of the specified parameters.

  procedure Secant_Single_Real_Predictor 
                ( x : in out Vector; prev_x : in Vector;
                  t : in out double_complex; prev_t,target : in double_complex;
                  h,tol : in double_float; pow : in positive := 1 );

  procedure Secant_Multiple_Real_Predictor 
                ( x : in out Solution_Array; prev_x : in Solution_Array;
                  t : in out double_complex; prev_t,target : in double_complex;
                  h,tol,dist_x : in double_float; pow : in positive := 1 );

  -- DESCRIPTION :
  --   Secant predictor for x and a linear predictor for t.

  -- ON ENTRY :
  --   x          the current approximation(s) for t;
  --   prev_x     the approximation(s) for a previous t;
  --   t          the current value of the continuation parameter;
  --   prev_t     is the previous value for t;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   pow        power of t in the homotopy towards 1.

  -- ON RETURN :
  --   x          the predicted approximation(s);
  --   t          the predicted value of the continuation parameter.

  procedure Secant_Single_Complex_Predictor 
                ( x : in out Vector; prev_x : in Vector;
                  t : in out double_complex; prev_t,target : in double_complex;
                  h,tol,dist_t : in double_float; trial : in natural );

  procedure Secant_Multiple_Complex_Predictor 
                ( x : in out Solution_Array; prev_x : in Solution_Array;
                  t : in out double_complex; prev_t,target : in double_complex; 
                  h,tol,dist_x,dist_t : in double_float; 
                  trial : in natural);

  -- DESCRIPTION :
  --   Secant predictor for x and complex predictor for t.

  -- ON ENTRY :
  --   x          the current approximation(s) for t;
  --   prev_x     the approximation(s) for a previous t;
  --   t          the current value of the continuation parameter;
  --   prev_t     is the previous value for t;
  --   target     is the target value for t;
  --   h          is the step size;
  --   tol        tolerance to decide when two numbers are equal;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   dist_t     t must keep a distance to the target;
  --   trial      indicates the number of trials for starting out of
  --              the previous value for t.

  -- ON RETURN :
  --   x          the predicted approximation(s);
  --   t          the predicted value of the continuation parameter.

  procedure Secant_Circular_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out double_complex; theta : in out double_float;
                  prev_t,t0_min_target,target : in double_complex;
                  h,tol : in double_float );

  -- DESCRIPTION :
  --   Secant predictor for x and circular predictor for t, around target.

  -- NOTE : This is the link between t and theta :
  --   t = target + t0_min_target * ( cos(theta) + i sin(theta) )

  -- ON ENTRY :
  --   x,prev_x,t,prev_t,target as before;
  --   theta      the angle for t.
  --   t0_min_target is t0-target, where t0 is the start point;
  --   h          is step size for theta !!

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter;
  --   theta      the predicted angle.

  procedure Secant_Geometric_Predictor
               ( x : in out Vector; prev_x : in Vector;
                 t : in out double_complex; prev_t,target : in double_complex;
                 h,tol : in double_float );

  -- DESCRIPTION :
  --   Secant predictor for x and a geometric predictor for t.

  -- ON ENTRY :
  --   x          the current approximation(s) for t;
  --   prev_x     the approximation(s) for a previous t;
  --   t          the current value of the continuation parameter;
  --   prev_t     is the previous value for t;
  --   target     is the target value for t;
  --   h          ratio between two consecutive distance to target, 0<h<1;
  --   tol        tolerance to decide when t = target;

  -- ON RETURN :
  --   x          the predicted approximation;
  --   t          the predicted value of the continuation parameter.

  generic

    with function Norm ( x : Vector ) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Single_Real_Predictor 
                ( x : in out Vector; t : in out double_complex;
                  target : in double_complex; h,tol : in double_float;
                  pow : in positive := 1 );

  generic

    with function Norm ( x : Vector ) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Multiple_Real_Predictor
                ( x : in out Solution_Array; t : in out double_complex;
                  target : in double_complex; h,tol,dist_x : in double_float;
                  nsys : in out natural; pow : in positive := 1 );

  -- DESCRIPTION :
  --   Tangent predictor for x and a linear predictor for t.

  -- ON ENTRY :
  --   x          current approximation for the solution;
  --   t          current value of the continuation parameter;
  --   target     target value for the continuation parameter;
  --   h          steplength;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > dist_x, for k in 1..n;
  --   nsys       must be initally equal to zero, used for counting;
  --   pow        power of t in the homotopy.

  -- ON RETURN :
  --   x          predicted approximation for the solution;
  --   t          new value of the continuation parameter;
  --   nsys       the number of linear systems solved.

  generic

    with function Norm ( x : Vector ) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Single_Complex_Predictor
                ( x : in out Vector; t : in out double_complex; 
                  target : in double_Complex;
                  h,tol,dist_t : in double_float; trial : in natural );

  generic

    with function Norm ( x : Vector) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Multiple_Complex_Predictor
                ( x : in out Solution_Array; t : in out double_complex;
                  target : in double_complex;
                  h,tol,dist_x,dist_t : in double_float;
	          trial : in natural; nsys : in out natural );

  -- DESCRIPTION :
  --   Tangent predictor for x and a complex predictor for t.

  -- ON ENTRY :
  --   x          current approximation for the solution;
  --   t          current value of the continuation parameter;
  --   target     target value for the continuation parameter;
  --   h          steplength;
  --   tol        tolerance to decide when two numbers are equal;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > d, for k in 1..n;
  --   dist_t     t must keep distance to the target;
  --   trial      indicates the number of trials for starting out of
  --              the previous value for t;
  --   nsys       must be initially equal to zero.

  -- ON RETURN :
  --   x          predicted approximation for the solution;
  --   t          new value of the continuation parameter;
  --   nsys       the number of linear systems solved.

  generic

    with function Norm ( x : Vector) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Circular_Predictor
                ( x : in out Vector; t : in out double_complex;
                  theta : in out double_float;
                  t0_min_target,target : in double_complex;
                  h,tol : in double_float );

  -- DESCRIPTION :
  --   This is a tangent predictor for x and a circular predictor for t
  --   For information on the parameters, see Secant_Circular_Predictor.

  generic

    with function Norm ( x : Vector) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Tangent_Geometric_Predictor
               ( x : in out Vector; t : in out double_complex;
                 target : in double_complex; h,tol : in double_float );

  -- DESCRIPTION :
  --   Tangent predictor for x and a geometric predictor for t.
  --   For information on the parameters, see Secant_Geometric_Predictor.

  generic

    with function Norm ( x : Vector) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Hermite_Single_Real_Predictor
                ( x : in out Vector; prev_x : in Vector;
                  t : in out double_complex; prev_t,target : in double_complex;
                  v : in out Vector; prev_v : in Vector;
                  h,tol : in double_float; pow : in positive := 1 );

  -- DESCRIPTION :
  --   Third-order extrapolation based on previous values of the solution
  --   paths along with corresponding first derivatives.

 -- ON ENTRY :
  --   x          current approximation for the solution;
  --   prev_x     previous approximation of the solution at prev_t;
  --   t          current value of the continuation parameter;
  --   prev_t     previous value of the continuation parameter;
  --   target     target value for the continuation parameter;
  --   v          will be used at work space;
  --   prev_v     direction of the path at prev_t;
  --   h          steplength;
  --   tol        tolerance to decide when t = target;
  --   dist_x     for all i /= j : |x(i)(k) - x(j)(k)| > dist_x, for k in 1..n;
  --   pow        power of t in the homotopy.

  -- ON RETURN :
  --   x          predicted approximation for the solution;
  --   t          new value of the continuation parameter;
  --   v          direction of the path computed for value of t on entry.

end Predictors;
