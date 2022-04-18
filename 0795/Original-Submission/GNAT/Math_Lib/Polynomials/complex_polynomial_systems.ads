with Vectors;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Vectors_of_Vectors;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package Complex_Polynomial_Systems is 

-- DESCRIPTION :
--   This package provides polynomial systems as vectors of polynomials,
--   added with data structures for evaluation and corresponding functions.

-- DATA STRUCTURES :

  package Vectors_of_Complex_Multivariate_Polynomials is
     new Vectors
       (Poly,Null_Poly,Clear,Copy,Equal,"+","-","-","*",
        Plus_Poly,Min_Poly,Min_Poly,Mult_Poly);

  type Poly_Sys is new Vectors_of_Complex_Multivariate_Polynomials.Vector;
  type Link_to_Poly_Sys is access Poly_Sys;

  type Eval_Poly_Sys is array ( integer range <> ) of Eval_Poly;
  type Eval_Coeff_Poly_Sys is array ( integer range <> ) of Eval_Coeff_Poly;

-- CREATORS :

  function Create ( p : Poly_Sys ) return Eval_Poly_Sys;
  function Create ( p : Poly_Sys ) return Eval_Coeff_Poly_Sys;

-- ADDITIONAL ARITHMETICAL OPERATIONS :

  function "*" ( a : double_complex; p : Poly_Sys ) return Poly_Sys;  
                                                            -- return a*p;
  function "*" ( p : Poly_Sys; a : double_complex ) return Poly_Sys;
                                                            -- return p*a;
  procedure Mult_Cmplx ( p : in out Poly_Sys; a : in double_complex );
                                                            -- p := a*p;

  function Eval ( p : Poly_Sys; x : double_complex; i : natural )
                return Poly_Sys;

  function Eval ( p : Poly_Sys;      x : Vector ) return Vector;
  function Eval ( p : Eval_Poly_Sys; x : Vector ) return Vector;
  function Eval ( p : Eval_Coeff_Poly_Sys;
                  c : Complex_Vectors_of_Vectors.Vector; x : Vector )
                return Vector;

  function  Diff ( p : Poly_Sys; i : natural ) return Poly_Sys;
  procedure Diff ( p : in out Poly_Sys; i : in natural );

-- DESTRUCTORS :

  procedure Clear ( p : in out Eval_Poly_Sys );
  procedure Clear ( p : in out Eval_Coeff_Poly_Sys );
  procedure Clear ( p : in out Link_to_Poly_Sys );

end Complex_Polynomial_Systems;
