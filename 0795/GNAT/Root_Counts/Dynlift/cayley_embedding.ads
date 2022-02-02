with Integer_Vectors;                    use Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Simplices,Triangulations;           use Simplices,Triangulations;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;

package Cayley_Embedding is

-- DESCRIPTION :
--   This procedure provides utilities to implement the embedding
--   of a tuple of supports to use the Cayley trick.
--   The novelty with this implementation is that the additional
--   coordinates to the points are added in front.

-- UTILITIES TO APPLY BEFORE THE CONSTRUCTION OF THE TRIANGULATION :

  function Embedding_Before_Lifting
               ( supports : Array_of_Lists ) return List;

  -- DESCRIPTION :
  --   The supports will be embedded in higher dimensional space
  --   by adding (r-1) extra coordinates.
  --   The supports will be placed on the vertices of a simplex.
  --   This is done before the lifting is defined on the supports.

-- UTILITIES TO APPLY AFTER THE CONSTRUCTION OF THE TRIANGULATION :

  function Extract ( vtp,n : natural; pts : List ) return List;

  -- DESCRIPTION :
  --   Extracts the points out of the list that are of the type
  --   indicated by vtp, i.e.: all points that belong to the list of 
  --   points placed on the vertex with number vtp.

  function Extract_Mixed_Cell
                ( n : natural; mix : Vector; s : Simplex ) return Mixed_Cell;

  -- DESCRIPTION :
  --   If the simplex determines a mixed cell, then this mixed cell
  --   will be returned.  Otherwise, the returned mixed cell will have
  --   an empty inner normal.

  function Extract_Mixed_Cells
                ( n : natural; mix : Vector; t : Triangulation )
                return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Given a triangulation of the large extended polytope,
  --   the mixed cells will be extracted.  The parameter
  --   n equals the dimension before the embedding and lifting.

  function Extract_non_Flat_Mixed_Cells
                ( n : natural; mix : Vector; t : Triangulation )
                return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Does the same as the routine just above, but stops when
  --   a simplex with normal (0,..,0,1) is encountered.

  procedure Deflate ( n : natural; l : in out List );

  -- DESCRIPTION :
  --   Removes the extra coordinates which have been introduced for
  --   the embedding.  The parameter n equals the dimension before
  --   the embedding and the lifting.

  procedure Deflate ( n : natural; mic : in out Mixed_Cell );

  -- DESCRIPTION :
  --   Removes the extra coordinates which have been introduced for
  --   the embedding.  The parameter n equals the dimension before
  --   the embedding and the lifting.

  procedure Deflate ( n : natural; mixsub : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Removes the extra coordinates which have been introduced for
  --   the embedding.  The parameter n equals the dimension before
  --   the embedding and the lifting.

end Cayley_Embedding;
