with Lists;
with Float_Vectors;                    use Float_Vectors;
with Float_Vectors_of_Vectors;

package Lists_of_Float_Vectors is

-- DESCRIPTION :
--   This package offers an abstraction for working with
--   lists of floating point vectors.

-- DATA STRUCTURE : a list of pointers to float vectors

  package Lists_of_Link_to_Float_Vectors is new Lists(Link_to_Vector);
  type List is new Lists_of_Link_to_Float_Vectors.List;

-- CONSTRUCTORS :

  function Deep_Create    ( v : Float_Vectors_of_Vectors.Vector ) return List;
  function Shallow_Create ( v : Float_Vectors_of_Vectors.Vector ) return List;

  function Deep_Create    ( l : List ) return Float_Vectors_of_Vectors.Vector;
  function Shallow_Create ( l : List ) return Float_Vectors_of_Vectors.Vector;

  -- DESCRIPTION :
  --   l := Create(v) equals v := Create(l).
  --   There is no sharing of pointers with a deep Create.
  --   With a shallow Create, both structure share the pointers.

  procedure Copy ( l1 : in List; l2 : in out List );

  -- DESCRIPTION :
  --   After Copy(l1,l2), Equal(l1,l2) holds.
  --   Of course, this is a deep copy, a shallow copy is given by l2 := l1.

  procedure Append ( first,last : in out List; v : in Vector );

  -- DESCRIPTION :
  --   The vector will be appended to the list first,
  --   last is a pointer to the last element of the list first.

  procedure Append_Diff ( first,last : in out List; v : in Vector );
  procedure Append_Diff ( first,last : in out List; v : in Link_to_Vector );

  -- DESCRIPTION :
  --   Only when v does not already belong to first, v will be added.

  procedure Deep_Concat    ( first,last : in out List; l : in List );
  procedure Shallow_Concat ( first,last : in out List; l : in List );

  -- DESCRIPTION :
  --   The list l will be concatenated to the list first,
  --   last is a pointer to the last element of the list first.
  --   With a deep concatenation, no pointers are shared.

  procedure Deep_Concat_Diff    ( first,last : in out List; l : in List );
  procedure Shallow_Concat_Diff ( first,last : in out List; l : in List );

  -- DESCRIPTION :
  --   Only those vectors of l will be concatenated that are not already
  --   in the list first.
  --   With a deep concatenation, no pointers are shared.

  procedure Remove ( l : in out List; x : in Vector );
  procedure Remove ( l : in out List; x : in Link_to_Vector );

  -- DESCRIPTION :
  --   Removes the point x from the list l.

  procedure Swap_to_Front ( l : in out List; x : in Vector );
  procedure Swap_to_Front ( l : in out List; x : in Link_to_Vector );

  -- DESCRIPTION :
  --   The point x belongs to the list l,
  --   Its content will be place in front of l and the first element
  --   of l will be swapped to the place of x.

-- SELECTORS :

  function Is_In ( l : List; v : Vector ) return boolean;
  function Is_In ( l : List; v : Link_to_Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector belongs to l.

  function Is_Equal ( l1,l2 : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both lists have the same vectors.

-- DESTRUCTOR :

  procedure Deep_Clear    ( l : in out List );
  procedure Shallow_Clear ( l : in out List );

  -- DESCRIPTION :
  --   Frees all allocated memory space.
  --   A deep clear deallocates also the points to the floating point vectors,
  --   while a shallow clear only removes the list structure.

end Lists_of_Float_Vectors;
