with text_io,Integer_Vectors;          use text_io,Integer_Vectors;
with Integer_Vectors_of_Vectors;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;
with Arrays_of_Integer_Vector_Lists;   use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;       use Integer_Mixed_Subdivisions;

package Mixed_Volume_Computation is

-- DESCRIPTION :
--   This package offers a number of routines for the computation
--   of the mixed volume of a system of polynomial equations.

-- UTILITIES :

  procedure Compute_Mixture ( supports : in out Array_of_Lists;
			      mix,perms : out Link_to_Vector ); 
  -- DESCRIPTION : 
  --   Computes the type of mixture of the supports of a system.

  -- ON ENTRY :
  --   supports   the supports of a polynomial system.

  -- ON RETURN :
  --   supports   a permuted array of supports, so that the same
  --              supports stand all toghether;
  --   mix        mix(k) indicates number of occurrencies of the kth support;
  --   perms      perms(k) gives the place of the kth support,
  --              after permutation to make supports correspond with mix.

  function Compute_Index ( k : natural; mix : Vector ) return natural;

  -- DESCRIPTION :
  --   Given k, an entry in the supports, the number this function returns
  --   indicates the number of different support, w.r.t. the type of  mixture.

  function Compute_Permutation
                  ( n : natural; mix : Vector; supports : Array_of_Lists )
                  return Link_to_Vector;

  -- DESCRIPTION :
  --   Given the type of mixture and the support, the permutation vector
  --   will be computed.

  -- ON RETURN :
  --   perms        perms(k) gives the place of the kth support,
  --                after permutation to make supports correspond with mix.

  function Typed_Lists ( mix : Vector; points : Array_of_Lists )
                       return Array_of_Lists;

  -- DESCRIPTION :
  --   Returns a tuple of lists where each list occurs only once,
  --   according to the given type of mixture.

  function Permute ( p : Poly_Sys; perm : Link_to_Vector ) return Poly_Sys;
  function Permute ( supports : Array_of_Lists ; perm : Link_to_Vector )
                   return Array_of_Lists;

  -- DESCRIPTION :
  --   Permutes the polynomials in the system or the supports,
  --   according to the vector perm.

-- MIXED VOLUME COMPUTATION, GIVEN A SUBDIVISION :

  function Mixed_Volume ( n : natural; mix : Vector;
                          mic : Mixed_Cell ) return natural;
  function Mixed_Volume ( n : natural; mix : Vector;
                          mixsub : Mixed_Subdivision ) return natural;

  -- DESCRIPTION :
  --   Computes the mixed volume based on a mixed cell and subdivision.
  --   When the cells are not fine enough, they will be refined but will
  --   be lost after returning the result.

  procedure Mixed_Volume ( n : in natural; mix : in Vector; 
                           mic : in out Mixed_Cell; mv : out natural );
  procedure Mixed_Volume ( n : in natural; mix : in Vector; 
                           mixsub : in out Mixed_Subdivision;
                           mv : out natural );

  -- DESCRIPTION :
  --   Computes the mixed volume based on a mixed cell and subdivision.
  --   When the cells are not fine enough, they will be refined by lifting.
  --   The refinement is stored in the subdivision field of the cells.

-- MIXED VOLUME COMPUTATIONS, GIVEN THE SUPPORTS :

  function Mixed_Volume ( n : natural; supports : Array_of_Lists )
			return natural;

  function Mixed_Volume ( file : file_type; n : natural;
                          supports : Array_of_Lists ) return natural;

  function Mixed_Volume ( n : natural; mix : Vector;
                          supports : Array_of_Lists ) return natural;

  function Mixed_Volume ( file : file_type; n : natural; mix : Vector;
                          supports : Array_of_Lists ) return natural;

  procedure Mixed_Volume ( n : in natural; mix : in Vector;
                           supports : in Array_of_Lists;
                           lifted : out Array_of_Lists;
                           mixsub : out Mixed_Subdivision; mv : out natural );

  procedure Mixed_Volume ( file : in file_type; n : in natural;
                           mix : in Vector; supports : in Array_of_Lists;
                           lifted : out Array_of_Lists;
                           mixsub : out Mixed_Subdivision; mv : out natural );

  -- DESCRIPTION :
  --   All these routines compute the mixed volume of support lists.

  -- ON ENTRY :
  --   file       if specified, then the mixed subdivision will be 
  --              written on file;
  --   n          the dimension of the system;
  --   mix        mix(k) is the number of times the kth support occurs;
  --   supports   the supports of a system of n polynomials in n unknowns.

  -- ON RETURN :
  --   lifted     array of listed points;
  --   mixsub     mixed subdivision used;
  --   mv         the mixed volume.

end Mixed_Volume_Computation;
