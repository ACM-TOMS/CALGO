with Vectors,Natural_Vectors;    use Natural_Vectors;

generic

  type coefftp is private;     -- type of the coefficients of the polynomials

  type Vector_of_coefftp is array ( integer range <> ) of coefftp;

  zero : coefftp;              -- the neutral element for "+"
                             
  with procedure clear ( a : in out coefftp );
  with procedure copy  ( a : in coefftp; b : in out coefftp );
  with function  equal ( a,b : coefftp ) return boolean;

  with function "+"     ( a,b : coefftp ) return coefftp;  -- return a+b;
  with function "-"     ( a : coefftp )   return coefftp;  -- return -a;
  with function "-"     ( a,b : coefftp ) return coefftp;  -- return a-b;
  with function "*"     ( a,b : coefftp ) return coefftp;  -- return a*b;
  with function convert ( n : natural )   return coefftp;  -- return coefftp(n);

  with procedure Plus_Coeff ( a : in out coefftp; b : in coefftp ); -- a := a+b;
  with procedure Min_Coeff  ( a : in out coefftp);                  -- a := -a;
  with procedure Min_Coeff  ( a : in out coefftp; b : in coefftp ); -- a := a-b;
  with procedure Mult_Coeff ( a : in out coefftp; b : in coefftp ); -- a := a*b;

  with function "<" ( v1,v2 : Link_to_Vector ) return boolean;
    -- return v1 < v2
  with function ">" ( v1,v2 : Link_to_Vector ) return boolean;
    -- return v1 > v2

package Multivariate_Polynomials is

-- DESCRIPTION :
--   This package provides a data abstraction for multivariate polynomials,
--   with coefficients over any ring, to be specified by instantiation.
--   Also the term ordening is determined by the generic parameters.
--   Two special data structures are provided for efficient evaluation.
--   Garbage collection is provided.  Therefore operations like clear,
--   equal and copy must be supplied when instantiating this package.

-- DATA STRUCTURES :

  type Degrees is new Natural_Vectors.Link_to_Vector;

  type Term is record
    cf : coefftp;    -- coefficient of the term
    dg : Degrees;    -- the degrees of the term
  end record;

  type Poly is private;
  type Eval_Poly is private;
  type Eval_Coeff_Poly is private;

  Null_Poly : constant Poly;

-- CONSTRUCTORS :

  function Create ( t : Term ) return Poly;
  function Create ( p : Poly ) return Eval_Poly;
  function Create ( p : Poly ) return Eval_Coeff_Poly;

  procedure Copy ( t1 : in Term; t2 : in out Term );    -- makes a deep copy
  procedure Copy ( p: in Poly; q : in out Poly );

-- SELECTORS :

  function Equal ( t1,t2 : Term )  return boolean;
  function Equal ( p,q : Poly )  return boolean;

  function Number_of_Unknowns ( p : Poly ) return natural;
  function Number_of_Terms    ( p : Poly ) return natural;

  function Degree ( p : Poly ) return integer;              -- return deg(p);
  function Degree ( p : Poly; i : integer ) return integer; -- return deg(p,xi);

  function "<" ( d1,d2 : Degrees ) return boolean;          -- return d1 < d2
  function ">" ( d1,d2 : Degrees ) return boolean;          -- return d1 > d2

  function Coeff ( p : Poly; d : Degrees ) return coefftp;
   -- Ex.: Coeff(c1*x^2+c2*x*y^3,(1 2))=c2;  Coeff(c1*x^2+c2,(1 0))=zero;
  function Coeff ( p : Poly ) return Vector_of_coefftp;
   -- returns the vector of coefficients, starting with highest degree

-- ARITHMETICAL OPERATIONS :

  function "+" ( p : Poly; t : Term ) return Poly;      -- return p+t;
  function "+" ( t : Term; p : Poly ) return Poly;      -- return t+p;
  function "+" ( p,q : Poly ) return Poly;              -- return p+q;
  function "-" ( p : Poly; t : Term ) return Poly;      -- return p-t;
  function "-" ( t : Term; p : Poly ) return Poly;      -- return t-p;
  function "-" ( p : Poly ) return Poly;                -- return -p;
  function "-" ( p,q : Poly ) return Poly;              -- return p-q;
  function "*" ( p : Poly; a : coefftp ) return Poly;   -- return a*p;
  function "*" ( a : coefftp; p : Poly ) return Poly;   -- return p*a;
  function "*" ( p : Poly; t : Term ) return Poly;      -- return p*t;
  function "*" ( t : Term; p : Poly ) return Poly;      -- return t*p;
  function "*" ( p,q : Poly ) return Poly;              -- return p*q;

  procedure Plus_Term ( p : in out Poly; t : in Term );    -- p := p + t;
  procedure Plus_Poly ( p : in out Poly; q : in Poly );    -- p := p + q;
  procedure Min_Term  ( p : in out Poly; t : in Term );    -- p := p - t;
  procedure Min_Poly  ( p : in out Poly );                 -- p := -p;
  procedure Min_Poly  ( p : in out Poly; q : in Poly );    -- p := p - q;
  procedure Mult_Coeff( p : in out Poly; a : in coefftp ); -- p := p * a;
  procedure Mult_Term ( p : in out Poly; t : in Term );    -- p := p * t;
  procedure Mult_Poly ( p : in out Poly; q : in Poly );    -- p := p * q;

  function Eval ( p : Poly; x : coefftp; i : integer ) return Poly;
     -- return p(x1,..,xi=x,..,xn);
     -- Number_of_Unknowns(Eval(p,x,i)) = Number_of_Unknowns(p)-1

  function Eval ( d : Degrees; c : coefftp; x : Vector_of_coefftp )
                return coefftp;  -- return c*x**d
  function Eval ( t : Term; c : coefftp; x : Vector_of_coefftp )
                return coefftp;  -- return c*x**d, with d = t.dg
  function Eval ( t : Term; x : Vector_of_coefftp ) return coefftp;

  function Eval ( p : Poly; x : Vector_of_coefftp ) return coefftp;
     -- return p(x);
  function Eval ( p : Poly; c,x : Vector_of_coefftp ) return coefftp;
     -- return p(c,x), with c = vector of coefficients for p

  function Eval ( p : Eval_Poly; x : Vector_of_coefftp ) return coefftp;
     -- return p(x);
  function Eval ( p : Eval_Coeff_Poly; c,x : Vector_of_coefftp ) return coefftp;
     -- return p(c,x), with c = vector of coefficients for p

  function  Diff ( p : Poly; i : integer ) return Poly; 
  procedure Diff ( p : in out Poly; i : in integer );
    -- symbolic differentiation w.r.t. the i-th unknown of p

  procedure Diff ( p : in Poly; i : in integer;
                   cp : out Eval_Coeff_Poly; m : out Vector_of_coefftp );
    -- evaluable coefficient polynomial of the partial derivative,
    -- with m the multiplication factors of the coefficients of p

-- ITERATORS : run through all terms of p and apply the generic procedure.

  generic
    with procedure process ( t : in out Term; continue : out boolean );
  procedure Changing_Iterator ( p : in out Poly );  -- t can be changed
  generic
    with procedure process ( t : in Term; continue : out boolean );
  procedure Visiting_Iterator ( p : in Poly );      -- t can only be read

-- DESTRUCTORS : deallocate memory.

  procedure Clear ( t : in out Term );
  procedure Clear ( p : in out Poly );
  procedure Clear ( p : in out Eval_Poly );
  procedure Clear ( p : in out Eval_Coeff_Poly );

private

  type Poly_Rep;
  type Eval_Poly_Rep;
  type Eval_Coeff_Poly_Rep;

  type Poly is access Poly_Rep;
  type Eval_Poly is access Eval_Poly_Rep;
  type Eval_Coeff_Poly is access Eval_Coeff_Poly_Rep;

  Null_Poly : constant Poly := null;

end Multivariate_Polynomials;
