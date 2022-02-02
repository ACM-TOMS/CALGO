with Floating_Point_Numbers;   use Floating_Point_Numbers;
with Complex_Vectors;          use Complex_Vectors;

package Complex_Norms is

-- DESCRIPTION :
--   This package offers some vector norms on complex vectors

  function Norm1 ( x : Vector ) return double_float;

  -- DESCRIPTION :
  --   returns the max |x(i)|, for i in x'range

  function Norm2 ( x : Vector ) return double_float;

  -- DESCRIPTION :
  --   returns sqrt(|sum {x(i)*x(i), for i in x'range}|)

 -- pragma inline(Norm1,Norm2);

end Complex_Norms;
