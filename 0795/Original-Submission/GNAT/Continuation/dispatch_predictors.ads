with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Matrices;                  use Complex_Matrices;
with Solutions,Continuation_Data;       use Solutions,Continuation_Data;

package Dispatch_Predictors is

-- DESCRIPTION :
--   This package provides generic predictors.

  generic

    with function Norm ( x : Vector ) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Single_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars;
                prev_x,prev_v : in Vector; v : in out Vector;
                prev_t,target : in double_complex;
                step,tol : in double_float; trial : in out natural );

  -- DESCRIPTION :
  --   Generic predictor for one solution.

  -- ON ENTRY :
  --   s        information about the current solution;
  --   p        parameters for the predictor;
  --   prev_x   previous solution component (only for secant);
  --   prev_t   previous value for t (only useful for secant);
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities;
  --   trial    number of consecutive trials (for complex predictor).

  -- ON RETURN :
  --   s        predicted value for solution.

  generic

    with function Norm ( x : Vector ) return double_float;
    with function dH ( x : Vector; t : double_complex ) return Vector;
          -- returns the derivatives of H(x,t) w.r.t. t in (x,t)
    with function dH ( x : Vector; t : double_complex ) return Matrix;
          -- returns the Jacobian matrix of H(x,t) at (x,t)

  procedure Multiple_Predictor
              ( s : in out Solu_Info_Array; p : in Pred_Pars;
                sa : in out Solution_Array; prev_sa : in Solution_Array; 
                t : in out double_complex; prev_t,target : in double_complex;
                step,tol,dist : in double_float; trial : in natural );

  -- DESCRIPTION :
  --   Generic predictor for an array of solutions.

  -- ON ENTRY :
  --   s        array with information of current solutions;
  --   sa       the current solutions;
  --   p        parameters for the predictor;
  --   prev_sa  previous solution component (only for secant);
  --   t        current value for continuation parameter;
  --   prev_t   previous value for t (only useful for secant);
  --   target   target value for continuation parameter;
  --   step     current step size;
  --   tol      tolerance for floating equalities;
  --   dist     tolerance for distance between solutions;
  --   trial    number of consecutive trials (for complex predictor).

  -- ON RETURN :
  --   sa       predicted values for solutions;
  --   t        predicted continuation parameter.

end Dispatch_Predictors;
