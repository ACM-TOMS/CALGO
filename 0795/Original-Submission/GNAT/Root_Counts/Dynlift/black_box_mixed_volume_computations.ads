with Integer_Vectors;                  use Integer_Vectors;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;
with Solutions;                        use Solutions;
with Arrays_of_Integer_Vector_Lists;   use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;       use Integer_Mixed_Subdivisions;

package Black_Box_Mixed_Volume_Computations is

  procedure Black_Box_Mixed_Volume_Computation
                 ( p : in Poly_Sys; mix : out Link_to_Vector;
                   lifsup : out Link_to_Array_of_Lists;
                   mixsub : out Mixed_Subdivision; mv : out natural );

  -- DESCRIPTION :
  --   Selects the appropriate algorithm to compute the mixed volume.

  -- ON ENTRY :
  --   p           polynomial system.

  -- ON RETURN :
  --   mix         type of mixture;
  --   lifsup      lifted supports of the system;
  --   mixsub      regular mixed-cell configuration;
  --   mv          mixed volume.

  procedure Black_Box_Polyhedral_Continuation
                 ( p : in Poly_Sys; mix : in Vector;
                   lifsup : in Array_of_Lists;
                   mixsub : in Mixed_Subdivision;
                   q : in out Poly_Sys; qsols : in out Solution_List );

  -- DESCRIPTION :
  --   Creates a random coefficient start system, based on the
  --   regular mixed-cell configuration.

  -- ON ENTRY :
  --   p           polynomial system;
  --   mix         type of mixture;
  --   lifsup      lifted supports of the system;
  --   mixsub      regular mixed-cell configuration;
  --   mv          mixed volume.

  -- ON RETURN :
  --   q           random coefficient start system;
  --   qsols       solutions of q.

end Black_Box_Mixed_Volume_Computations;
