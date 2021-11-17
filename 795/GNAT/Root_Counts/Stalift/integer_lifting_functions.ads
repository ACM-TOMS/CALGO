with Integer_Vectors;                   use Integer_Vectors;
with Integer_Vectors_of_Vectors;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;

package Integer_Lifting_Functions is

-- DESCRIPTION :
--   This package offers three different types of integer-valued
--   lifting functions: linear, polynomial, and point-wise.

-- LINEAR LIFTING :

  function Linear_Lift ( lf,v : Vector ) return Vector;
  function Linear_Lift ( lf : Vector; l : List ) return List;
  function Linear_Lift ( lf : Integer_Vectors_of_Vectors.Vector;
                         l : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   The last entry of the enlarged vector on return equals <lf,v>,
  --   for every vector v in the lists.

  -- REQUIRED : lf'range = v'range.

-- POLYNOMIAL LIFTING :

  function Polynomial_Lift ( lf : Poly; v : Vector ) return Vector;
  function Polynomial_Lift ( lf : Eval_Poly; v : Vector ) return Vector;
  function Polynomial_Lift ( lf : Poly; l : List ) return List;
  function Polynomial_Lift ( lf : Eval_Poly; l : List ) return List;
  function Polynomial_Lift ( lf : Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists;
  function Polynomial_Lift ( lf : Eval_Poly_Sys; l : Array_of_Lists )
                           return Array_of_Lists;

  -- DESCRIPTION :
  --   As lf is complex valued, the rounded real part is used as lifting.

-- RANDOM LIFTING :

  function Random_Lift ( lflow,lfupp : integer; v : Vector ) return Vector;
  function Random_Lift ( lflow,lfupp : integer; l : List ) return List;
  function Random_Lift ( lflow,lfupp : Vector; l : Array_of_Lists )
                       return Array_of_Lists;

  -- DESCRIPTION :
  --   The lifting values are random integers between lflow and lfupp.

-- RANDOM LINEAR LIFTING :

  function Random_Linear_Lift ( lflow,lfupp : integer; v : Vector ) 
                              return Vector;
  function Random_Linear_Lift ( lflow,lfupp : integer; l : List ) return List;
  function Random_Linear_Lift ( lflow,lfupp : Vector; l : Array_of_Lists )
                              return Array_of_Lists;

  -- DESCRIPTION :
  --   Linear lifting with random vectors of integers between lflow and lfupp.

-- POINT-WISE LIFTING :

  function Point_Lift ( lf : integer; v : Vector ) return Vector;
  function Point_Lift ( lf : Vector; l : List ) return List;
  function Point_Lift ( lf : Integer_Vectors_of_Vectors.Vector;
                        l : Array_of_Lists ) return Array_of_Lists;

  -- DESCRIPTION :
  --   The enlarged vector on return has as last component the lifting lf(k),
  --   where k is the position of the vector v in the lists.

  -- REQUIRED : Length_Of(l) = lf'length.

end Integer_Lifting_Functions;
