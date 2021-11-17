with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Float_Vectors;                      use Float_Vectors;
with Lists_of_Float_Vectors;             use Lists_of_Float_Vectors;
with Arrays_of_Float_Vector_Lists;       use Arrays_of_Float_Vector_Lists;
with Float_Faces_of_Polytope;            use Float_Faces_of_Polytope;
with Complex_Multivariate_Polynomials;   use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;

package Float_Lifting_Functions is

-- DESCRIPTION :
--   This package provides a suite of floating-point lifting functions.

-- RANDOM FLOATING-POINT LIFTING :

  function Random_Lift ( lflow,lfupp : double_float ) return double_float;
  function Random_Lift ( v : Vector; lflow,lfupp : double_float ) return Vector;
  function Random_Lift ( l : List; lflow,lfupp : double_float ) return List;
  function Random_Lift ( l : Array_of_Lists; lflow,lfupp : Vector )
                       return Array_of_Lists;
  -- DESCRIPTION :
  --   Random lifting values between lflow and lfupp.

-- LINEAR LIFTING FUNCTIONS :

  function Linear_Lift ( x,v : Vector ) return Vector;
  function Linear_Lift ( f : Face; v : Vector ) return Face;
  function Linear_Lift ( l : List; v : Vector ) return List;
  function Linear_Lift ( f : Faces; v : Vector ) return Faces;

  -- DESCRIPTION :
  --   Returns a linearly lifted vector, list or faces.

-- RANDOM FLOATING-POINT LINEAR LIFTING FUNCTIONS :

  function Random ( n : natural; lflow,lfupp : double_float ) return Vector;

  -- DESCRIPTION :
  --   Returns a vector of range 1..n with randomly generated numbers,
  --   in [lflow,lfupp].  Random linear lifting functions are provided
  --   by using this randomly generated vector.

-- POLYNOMIAL LIFTING FUNCTIONS :

  function Polynomial_Lift ( lf : Poly; x : Vector ) return Vector;
  function Polynomial_Lift ( lf : Eval_Poly; x : Vector ) return Vector;
  function Polynomial_Lift ( lf : Poly; l : List ) return List;
  function Polynomial_Lift ( lf : Eval_Poly; l : List ) return List;
  function Polynomial_Lift ( lf : Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists;
  function Polynomial_Lift ( lf : Eval_Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists;

end Float_Lifting_Functions;
