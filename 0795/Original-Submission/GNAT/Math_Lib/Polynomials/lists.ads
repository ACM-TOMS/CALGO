generic

  type Item is private;

package Lists is

-- DESCRIPTION :
--   This generic package allows to implement lists of items.

  type List is private;

  Null_List : constant List;

  Overflow     : exception;
  List_Is_Null : exception;

-- CONSTRUCTORS :

  procedure Construct ( i : in Item; l : in out List );

  -- DESCRIPTION :
  --   Adds the item i to the front of the list l.

  procedure Set_Head ( l : in out List; i : in Item);

  -- DESCRIPTION :
  --   Sets the first element in the list to item i.

  -- REQUIRED : not Is_Null(l).

  procedure Swap_Tail ( l1,l2 : in out List );

  -- DESCRIPTION :
  --   Swaps the tail of list l1 with the list l2.

  procedure Append ( first,last : in out List; i : in Item );

  -- DESCRIPTION :
  --   Appends the item i to the list, where first points to the first
  --   element and last to its last element.

  procedure Concat ( first,last : in out List; l : in List );

  -- DESCRIPTION :
  --   Concatenates the list l to the list first, where last points to
  --   the last element of the list.

  procedure Copy ( l1 : in List; l2 : in out List );

  -- DESCRIPTION :
  --   Makes a copy from the list l1 to the list l2.

-- SELECTORS :

  function Is_Equal ( l1,l2 : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both lists are equal.

  function Length_Of ( l : List ) return natural;

  -- DESCRIPTION :
  --   Returns the length of the list, i.e.: the number of elements.

  function Is_Null ( l : List ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the list is empty, false otherwise.

  function Head_Of ( l : List) return Item;

  -- DESCRIPTION :
  --   Returns the first element in the list.

  -- REQUIRED : not Is_Null(l).

  function Tail_Of ( l : List ) return List;

  -- DESCRIPTION :
  --   Returns the tail of the list l.

  -- REQUIRED : not Is_Null(l).

-- DESTRUCTOR :

  procedure Clear ( l : in out List );

  -- DESCRIPTION :
  --   Deallocates the memory occupied by the elements in the list.
    
private

  type Node;
  type List is access Node;
  Null_List : constant List := null;

end Lists;
