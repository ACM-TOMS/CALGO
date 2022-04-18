with Integer_Vectors,Float_Vectors;    use Integer_Vectors;
with Integer_Vectors_of_Vectors; 
with Arrays_of_Integer_Vector_Lists;   use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;       use Integer_Mixed_Subdivisions;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;

package Mixed_Coherent_Subdivisions is

-- DESCRIPTION :
--   A number of routines for constructing a mixed coherent subdivision
--   are offered, each with a different type of lifting function.

-- a polynomial system as lifting function :

  function Mixed_Coherent_Subdivision
               ( n : natural; mix : Vector; points : Array_of_Lists;
                 lift : Poly_Sys ) return Mixed_Subdivision;

  procedure Mixed_Coherent_Subdivision
               ( n : in natural; mix : in Vector; points : in Array_of_Lists;
                 lift : in Poly_Sys; lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

-- a user-defined lifting function :

  function Mixed_Coherent_Subdivision
               ( n : natural; mix : Vector; points : Array_of_Lists; 
                 linear : boolean; lift : Integer_Vectors_of_Vectors.Vector )
               return Mixed_Subdivision;

  procedure Mixed_Coherent_Subdivision
               ( n : in natural; mix : in Vector; points : in Array_of_Lists;
                 linear : in boolean;
                 lift : in Integer_Vectors_of_Vectors.Vector;
                 lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

-- a randomly generated lifting function :

  function Mixed_Coherent_Subdivision
               ( n : natural; mix : Vector; points : Array_of_Lists;
                 linear : boolean; low,upp : Vector )
               return Mixed_Subdivision;

  procedure Mixed_Coherent_Subdivision
	       ( n : in natural; mix : in Vector; points : in Array_of_Lists;
                 linear : in boolean; low,upp : in Vector;
                 lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Given a set of points and a lifting function,
  --   a subdivision of the polytope will be computed.

  -- ON ENTRY :
  --   n         the dimension of the vector space;
  --   mix       mix(k) indicates how many times the kth point set occurs;
  --   points    an array of all different point sets;
  --   linear    indicates wether a linear lifting should be used;
  --   lift      an array of lifting polynomials or an m dimensional
  --             array of vectors, where the length of the kth vector
  --             must equal the length of the kth support,
  --             when nonlinear, otherwise the length equals n.
  --   low,upp   lower and upper bounds for random lifting.

  -- ON RETURN :
  --   lifted    the lifted points which can later be used for lifting
  --             the polynomial system;
  --   nbsucc    the number of successful face-face combinations that
  --             have been computed;
  --   nbfail    the number of unsuccessful face-face combinations;
  --   mixsub    the mixed subdivision of the polytope, defined as
  --             lower hull of the lifted points.

end Mixed_Coherent_Subdivisions;
