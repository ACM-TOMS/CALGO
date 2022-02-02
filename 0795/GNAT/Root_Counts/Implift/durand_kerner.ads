with Floating_Point_Numbers;   use Floating_Point_Numbers;
with Complex_Vectors;          use Complex_Vectors;

generic 

  with procedure Write ( step : in natural; z,res : in Vector );

  -- DESCRIPTION :
  --   This routine allows to write intermediate results after each iteration,
  --   such as the step number, the approximations z and the residuals res.
  --   If no output is wanted, supply an empty body for Write.

procedure Durand_Kerner ( p : in Vector; z,res : in out Vector;
                          maxsteps : in natural; eps : in double_float;
                          nb : out natural );

-- DESCRIPTION :
--   This routine computes all roots of a given polynomial
--   in one unknown, applying the method of Durand Kerner.

-- ON ENTRY :
--   p        the polynomial defined by
--             p[k] + p[k+1]*x + p[k+2]*x^2 + .. + p[k+n]*x^n,
--            with k = p'first;
--   z        initial approximations for the roots;
--   res      the residuals of the roots;
--   maxsteps is the maximum number of steps that are allowed;
--   eps      the required accuracy

-- ON RETURN :
--   z        the computed roots;
--   res      the residuals of the roots;
--   nb       the number of steps.
