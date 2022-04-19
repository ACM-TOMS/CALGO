with Vectors;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Vectors_of_Vectors;
with Complex_Multivariate_Laurent_Polynomials;
 use Complex_Multivariate_Laurent_Polynomials;

package Complex_Laurent_Polynomial_Systems is 

-- DESCRIPTION :
--   This package provides polynomial systems as vectors of polynomials,
--   added with data structures for evaluation and corresponding functions.

-- DATA STRUCTURES :

  package Vectors_of_Complex_Multivariate_Laurent_Polynomials is
     new Vectors
       (Poly,Null_Poly,Clear,Copy,Equal,"+","-","-","*",
        Plus_Poly,Min_Poly,Min_Poly,Mult_Poly);

  type Laur_Sys is
         new Vectors_of_Complex_Multivariate_Laurent_Polynomials.Vector;
  type Link_to_Laur_Sys is access Laur_Sys;

  type Eval_Laur_Sys is array ( integer range <> ) of Eval_Poly;
  type Eval_Coeff_Laur_Sys is array ( integer range <> ) of Eval_Coeff_Poly;

-- CREATORS :

  function Create ( p : Laur_Sys ) return Eval_Laur_Sys;
  function Create ( p : Laur_Sys ) return Eval_Coeff_Laur_Sys;

-- ADDITIONAL ARITHMETIC OPERATIONS :

  function "*" ( a : double_complex; p : Laur_Sys ) return Laur_Sys;
                                                             -- return a*p;
  function "*" ( p : Laur_Sys; a : double_complex ) return Laur_Sys;
                                                             -- return p*a;
  procedure Mult_Cmplx ( p : in out Laur_Sys; a : in double_complex );
                                                             -- p := a*p;

  function Eval ( p : Laur_Sys; x : double_complex; i : natural )
                return Laur_Sys;

  function Eval ( p : Laur_Sys;      x : Vector ) return Vector;
  function Eval ( p : Eval_Laur_Sys; x : Vector ) return Vector;
  function Eval ( p : Eval_Coeff_Laur_Sys;
                  c : Complex_Vectors_of_Vectors.Vector; x : Vector )
                return Vector;

  function  Diff ( p : Laur_Sys; i : natural ) return Laur_Sys;
  procedure Diff ( p : in out Laur_Sys; i : in natural );

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Laur_Sys );
  procedure Clear ( p : in out Eval_Coeff_Laur_Sys );

end Complex_Laurent_Polynomial_Systems;
