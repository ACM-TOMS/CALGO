with text_io,Integer_Vectors;             use text_io,Integer_Vectors;
with Lists_of_Integer_Vectors;            use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;      use Arrays_of_Integer_Vector_Lists;
with Complex_Polynomial_Systems;          use Complex_Polynomial_Systems;

package Drivers_for_Vertex_Points is

-- DESCRIPTION :
--   This package provides two drivers for extracting the vertex
--   point out of a tuple of point lists.

  procedure Vertex_Points
                ( file : in file_type; l : in out List );
  procedure Vertex_Points 
                ( file : in file_type; l : in out Array_of_Lists );
  procedure Vertex_Points 
                ( file : in file_type; mix : in Link_to_Vector;
                  l : in out Array_of_Lists );

  -- DESCRIPTION :
  --   Reduces the lists to the lists of vertex points.

  -- REQUIRED :
  --   If the type of mixture (mix) is provided, then the tuple of lists
  --   must be sorted according to this vector mix.

  -- ON ENTRY :
  --   file       for writing diagnostics and statistics;
  --   mix        number of different lists in the tuple l,
  --              if not provided, then it will be assumed that all lists
  --              are different from each other;
  --   l          (tuple of) list(s).

  -- ON RETURN :
  --   l          (tuple of) list(s) with nothing but vertex points.

  procedure Vertex_Points
                ( file : in file_type; p : in out Poly_Sys );
  procedure Vertex_Points
                ( file : in file_type; mix : in Link_to_Vector;
                  p : in out Poly_Sys );

  -- DESCRIPTION :
  --   Reduces the supports of the polynomials to their vertex points.
  --   Merely a driver to the procedures listed above.

end Drivers_for_Vertex_Points;
