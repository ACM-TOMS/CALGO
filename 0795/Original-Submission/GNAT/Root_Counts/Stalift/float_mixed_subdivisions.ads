with Lists;
with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Float_Vectors;                    use Float_Vectors;
with Lists_of_Float_Vectors;           use Lists_of_Float_Vectors;
with Arrays_of_Float_Vector_Lists;     use Arrays_of_Float_Vector_Lists;

package Float_Mixed_Subdivisions is

-- DESCRIPTION :
--   This package enables working with regular mixed subdivisions.

-- DATA STRUCTURES :

  type Mixed_Subdivision;              -- list of mixed cells
  type Link_to_Mixed_Subdivision is access Mixed_Subdivision;

  type Mixed_Cell is record
    nor : Link_to_Vector;              -- inner normal to the facet
    pts : Link_to_Array_of_Lists;      -- points that span the cell
    sub : Link_to_Mixed_Subdivision;   -- subdivision of the cell
  end record;

  package Lists_of_Mixed_Cells is new Lists(Mixed_Cell);
  type Mixed_Subdivision is new Lists_of_Mixed_Cells.List;

-- CREATORS :

  function Create ( pts : Array_of_Lists; nor : Vector; tol : double_float )
                  return Mixed_Cell;
  function Create ( pts : Array_of_Lists; nors : List; tol : double_float )
                  return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Creates the mixed cell(s) of those points whose inner product
  --   with the given normal(s) is minimal.
  --   The parameter tol is the tolerance on the precision.

  function Create ( pts : Array_of_Lists; mixsub : Mixed_Subdivision;
                    tol : double_float ) return Mixed_Subdivision;

  -- DESCRIPTION :
  --   Takes the normals of the cells in the given mixed subdivision
  --   and creates the mixed cells by selecting the points whose inner
  --   product with the normals are minimal.

  procedure Update ( pts : in Array_of_Lists; nor : in Vector;
                     mixsub,mixsub_last : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Given a tuple of points and a normal, the mixed subdivision will
  --   be updated: either an existing cell will get additional points,
  --   if the normal already occurs in the subdivision, or otherwise,
  --   a new cell will be created and appended to the mixed subdivision.

-- CONSTRUCTORS :

  procedure Copy ( mic1 : in Mixed_Cell; mic2 : in out Mixed_Cell );
  procedure Copy ( mixsub1 : in Mixed_Subdivision; 
                   mixsub2 : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Makes a deep copy of the cells and the subdivisions.

  procedure Append_Diff ( first,last : in out Mixed_Subdivision;
                          mic : in Mixed_Cell );

  -- DESCRIPTION :
  --   Appends a mixed cell to the list of cells first, where
  --   last points to the last element of the list first.
  --   The suffix _Diff means that only when the cell does not already
  --   belong to the list first, it will be appended.

  procedure Concat_Diff ( first,last : in out Mixed_Subdivision;
                          mixsub : in Mixed_Subdivision );

  -- DESCRIPTION :
  --   Concatenates all cells in mixsub to the list of cells first,
  --   last is a pointer to the last cell in first.
  --   The suffix _Diff means that only when those cells that do not already
  --   belong to the list first will be appended.

  procedure Construct ( mixsub : in Mixed_Subdivision;
                        first : in out Mixed_Subdivision );
  procedure Construct_Diff ( mixsub : in Mixed_Subdivision;
                             first : in out Mixed_Subdivision );

  -- DESCRIPTION :
  --   Constructs all cells in the mixed subdivision to the front of
  --   the list first.
  --   The suffix _Diff means that only when those cells that do not already
  --   belong to the list first will be constructed to first.

-- SELECTORS :

  function Is_Equal ( mic1,mic2 : Mixed_Cell ) return boolean;
  function Is_Equal ( mixsub1,mixsub2 : Mixed_Subdivision ) return boolean;
  function Is_Equal ( mixsub1,mixsub2 : Link_to_Mixed_Subdivision )
                    return boolean;

  -- DESCRIPTION :
  --   Returns true when two mixed cells and mixed subdivisions are equal.

  function Is_In ( mixsub : Mixed_Subdivision; normal : Vector ) return boolean;
  function Is_In ( mixsub : Mixed_Subdivision; mic : Mixed_Cell )
                 return boolean;

  -- DESCRIPTION :
  --   Returns true if normal or cell belongs to the mixed subdivision,
  --   otherwise false is returned.  When the whole mixed cell is given,
  --   then not only the normal, but also the points will be checked.

-- DESTRUCTORS :

  procedure Deep_Clear ( mic : in out Mixed_Cell );
  procedure Deep_Clear ( mixsub : in out Mixed_Subdivision );
  procedure Deep_Clear ( mixsub : in out Link_to_Mixed_Subdivision );
  procedure Shallow_Clear ( mic : in out Mixed_Cell );
  procedure Shallow_Clear ( mixsub : in out Mixed_Subdivision );
  procedure Shallow_Clear ( mixsub : in out Link_to_Mixed_Subdivision );

  -- DESCRIPTION :
  --   The memory space allocated will be freed.
  --   A shallow clear only destroys the list structures,
  --   while with a deep clear, also the contents of the lists are freed.

end Float_Mixed_Subdivisions;
