with text_io;                          use text_io;
with Integer_Vectors,Triangulations;   use Integer_Vectors,Triangulations;
with Integer_Mixed_Subdivisions;       use Integer_Mixed_Subdivisions;

procedure Driver_for_Minkowski_Polynomials
                ( file : in file_type;
                  n : in natural; mix : in Vector; t : in Triangulation;
                  alltri : in boolean; mixsub : out Mixed_Subdivision );

-- DESCRIPTION :
--   Driver for the computation of the Minkowski-polynomial.

-- ON ENTRY :
--   file         to write all results on;
--   n            dimension before lifting and embedding;
--   mix          type of mixture;
--   t            triangulation of the Cayley polytope;
--   alltri       true when all triangulations are wanted, false otherwise.

-- ON OUTPUT :
--   mixed        mixed subdivision, corresponding the type of mixture.
