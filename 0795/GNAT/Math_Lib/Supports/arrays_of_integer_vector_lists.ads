with Lists_of_Integer_Vectors;      use Lists_of_Integer_Vectors;

package Arrays_of_Integer_Vector_Lists is

-- DESCRIPTION :
--   This package offers an abstraction for working with
--   arrays of lists of integer vectors.

-- DATA STRUCTURES :

  type Array_of_Lists is array ( integer range <> ) of List;
  type Link_to_Array_of_Lists is access Array_of_Lists;

-- CONSTRUCTOR :

  procedure Copy ( l1 : in Array_of_Lists; l2 : in out Array_of_Lists );

  -- DESCRIPTION :
  --   After Copy(l1,l2), Equal(l1,l2) holds.
  --   Of course, this is a deep copy, in constrast to l2 := l1.

-- SELECTORS :

  function Is_Equal ( l1,l2 : Array_of_Lists ) return boolean;
  function Is_Equal ( l1,l2 : Link_to_Array_of_Lists ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both arrays have the same lists.

  function Length_Of ( l : Array_of_Lists ) return natural;

  -- DESCRIPTION :
  --   Returns the sum of all lengths of the lists in l.

-- DESTRUCTORS :

  procedure Deep_Clear    ( l : in out Array_of_Lists );
  procedure Shallow_Clear ( l : in out Array_of_Lists );
  procedure Deep_Clear    ( l : in out Link_to_Array_of_Lists );
  procedure Shallow_Clear ( l : in out Link_to_Array_of_Lists );

  -- DESCRIPTION :
  --   Frees allocated memory space.
  --   With a deep clear, also the content of the lists are cleared,
  --   while with a shallow clear, only the lists structures will be
  --   destroyed, the points in the lists will remain.

end Arrays_of_Integer_Vector_Lists;
