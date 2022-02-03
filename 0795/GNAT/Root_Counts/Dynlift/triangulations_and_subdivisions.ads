with Integer_Vectors;               use Integer_Vectors;
with Simplices,Triangulations;      use Simplices,Triangulations;
with Integer_Mixed_Subdivisions;    use Integer_Mixed_Subdivisions;

package Triangulations_and_Subdivisions is

-- DESCRIPTION :
--   This package offers conversion operations from triangulations
--   into mixed subdivisions and vice versa.

  function Deep_Create    ( n : natural; s : Simplex ) return Mixed_Cell;
  function Shallow_Create ( n : natural; s : Simplex ) return Mixed_Cell;

  -- DESCRIPTION :
  --   Creates a mixed cell from the given simplex.  A deep create makes
  --   a copy of the points in the simplex, whereas a shallow create only
  --   copies the points to the points in the simplex.

  function Deep_Create    ( n : natural; t : Triangulation )
                          return Mixed_Subdivision;
  function Shallow_Create ( n : natural; t : Triangulation )
                          return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Creates a mixed subdivision from a given triangulation.

  function Deep_Create    ( n : natural; flatnor : Vector; t : Triangulation )
                          return Mixed_Subdivision;
  function Shallow_Create ( n : natural; flatnor : Vector; t : Triangulation )
                          return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Stops the conversion as soon as the flattening normal is encountered.
  --   Cells with the same inner normal are merged.

  function Non_Flat_Deep_Create    ( n : natural; t : Triangulation )
                                   return Mixed_Subdivision;
  function Non_Flat_Shallow_Create ( n : natural; t : Triangulation )
                                   return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Converts only the cells in the triangulation that are not flat.

  function Deep_Create    ( n : natural; mixsub : Mixed_Subdivision )
                          return Triangulation;
  function Shallow_Create ( n : natural; mixsub : Mixed_Subdivision )
                          return Triangulation;

  -- DESCRIPTION :
  --   Creates a triangulation from the mixed subdivision.

  -- REQUIRED :
  --   The subdivision must be a triangulation!

end Triangulations_and_Subdivisions;
