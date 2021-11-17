with Complex_Vectors,Complex_Matrices;  use Complex_Vectors,Complex_Matrices;
with Complex_Vectors_of_Vectors;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;

package Jacobi_Matrices is

-- DESCRIPTION:
--   This package contains routines for how to deal with
--   Jacobi Matrices of Complex Polynomial Systems.

-- DATA STRUCTURES :

  type Jacobi      is array ( integer range <>, integer range <> ) of Poly;
  type Eval_Jacobi is array ( integer range <>, integer range <> ) of Eval_Poly;

  type Eval_Coeff_Jacobi
             is array ( integer range <>, integer range <> ) of Eval_Coeff_Poly;
  type Mult_Factors 
             is array ( integer range <>, integer range <> ) of Link_to_Vector;
 
  -- USAGE :
  --   p : Poly_Sys;
  --   j : Jacobi := Create(p);
  --  =>  j(i,j) = Diff(p(i),j)

-- CREATORS :

  function Create ( p : Poly_Sys ) return Jacobi;

  -- REQUIRED :
  --   The number of the unknowns of each polynomial must be the same

  function Create ( j : Jacobi ) return Eval_Jacobi;

  procedure Create ( p : Poly_Sys; 
                     j : out Eval_Coeff_Jacobi; m : out Mult_Factors );

-- OPERATIONS :

  function  Equal ( j1,j2 : Jacobi ) return boolean;
  procedure Copy  ( j1 : in Jacobi; j2 : in out Jacobi );

  function "+" ( j1,j2 : Jacobi )          return Jacobi;      -- return j1+j2;
  function "-" ( j : Jacobi )              return Jacobi;      -- return -j;
  function "-" ( j1,j2 : Jacobi )          return Jacobi;      -- return j1-j2;

  procedure Plus_Jacobi ( j1 : in out Jacobi; j2 : in Jacobi ); -- j1 := j1+j2;
  procedure Min_Jacobi  ( j  : in out Jacobi );                 -- j  := -j;
  procedure Min_Jacobi  ( j1 : in out Jacobi; j2 : in Jacobi ); -- j1 := j1-j2;

  function Eval ( j : Jacobi;      x : Vector ) return matrix; -- return j(x);
  function Eval ( j : Eval_Jacobi; x : Vector ) return matrix; -- return j(x);

  function Eval ( j : Eval_Coeff_Jacobi; m : Mult_Factors; 
                  c : Complex_Vectors_of_Vectors.Vector; x : Vector )
                return Matrix;

    -- returns j(c,x) with c the coefficients of the original polynomials
  
-- DESTRUCTORS :

  procedure Clear ( j : in out Jacobi );
  procedure Clear ( j : in out Eval_Jacobi );
  procedure Clear ( j : in out Eval_Coeff_Jacobi );

  procedure Clear ( m : in out Mult_Factors );

end Jacobi_Matrices;
