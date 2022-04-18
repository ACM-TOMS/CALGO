with Symmetry_Group;                    use Symmetry_Group;
with Integer_Vectors;                   use Integer_Vectors;
with Integer_Mixed_Subdivisions;
with Float_Mixed_Subdivisions;

package Generating_Mixed_Cells is

-- DESCRIPTION :
--   Given a symmetric subdivision, the following routines return the
--   generating cells in the subdivision.

  function Generating_Cells
              ( v,w : List_of_Permutations; mix : Vector;
                mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision )
              return Integer_Mixed_Subdivisions.Mixed_Subdivision;

  -- DESCRIPTION :
  --   Extracts the generating cells from the (G,V,W)-symmetric subdivision,

  function Generating_Cells
              ( v,w : List_of_Permutations; mix : Vector;
                mixsub : Float_Mixed_Subdivisions.Mixed_Subdivision )
              return Float_Mixed_Subdivisions.Mixed_Subdivision;

  -- DESCRIPTION :
  --   Extracts the generating cells from the (G,V,W)-symmetric subdivision,

  function Generating_Cells
              ( mixsub : Integer_Mixed_Subdivisions.Mixed_Subdivision )
              return Integer_Mixed_Subdivisions.Mixed_Subdivision;

  -- DESCRIPTION :
  --   Extracts the generating mixed cells from mixsub, under the
  --   assumption that it is invariant under the full permutation group.

  function Generating_Cells
              ( mixsub : Float_Mixed_Subdivisions.Mixed_Subdivision )
              return Float_Mixed_Subdivisions.Mixed_Subdivision;

  -- DESCRIPTION :
  --   Extracts the generating mixed cells from mixsub, under the
  --   assumption that it is invariant under the full permutation group.

end Generating_Mixed_Cells;
