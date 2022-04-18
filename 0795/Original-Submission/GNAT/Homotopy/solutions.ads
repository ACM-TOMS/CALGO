with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Lists;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;

package Solutions is

-- DESCRIPTION :
--   This package provides an abstraction of a list and an array
--   of solutions.

-- DATA STRUCTURES :
 
  type Solution ( n : natural ) is record
    t : double_complex;         -- continuation parameter t
    m : natural;                -- multiplicity of the solution
    v : Vector(1..n);           -- the solution
    err : double_float;         -- error = |correction| from Newton
    rco : double_float;         -- inverse of condition number
    res : double_float;         -- norm of residual vector
  end record;

  type Link_to_Solution is access Solution;

  package List_of_Solutions is new Lists (Item => Link_to_Solution);
  type Solution_List is new List_of_Solutions.List;

  type Solution_Array is array ( positive range <> ) of Link_to_Solution;

-- CREATORS :

  function Create ( sl : Solution_List ) return Solution_Array;
  function Create ( sa : Solution_Array ) return Solution_List;

  -- DESCRIPTION :
  --   Allows the transition from a list to an array and vice versa.

-- SELECTORS :

  function Equal ( s1,s2 : Solution; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if for each component 
  --   |s1.v(i)-s2.v(i)|/|s1.v(i)| < tol, for i in s1.v'range.

  function Equal ( s1,s2 : Solution_List;  tol : double_float ) return boolean;
  function Equal ( s1,s2 : Solution_Array; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if both lists of arrays are equal to each other, upon
  --   the given tolerance for the relative error.

  procedure Equals ( sols : in out Solution_List; flag : in natural;
                     tol : in double_float; same : out boolean );
  -- DESCRIPTION :
  --   The solutions that are equal to each other are marked with a flag.

  procedure Equals ( sa : in Solution_Array; x : in Vector; i : in natural;
                     tol : in double_float; j : in out natural );
  -- DESCRIPTION :
  --   Compares the first i-1 vectors in sa with x.

  -- ON ENTRY :
  --   sa        a solution array, containing at least i-1 elements;
  --   x         a vector;
  --   i         an index, normally the entry of x in sa;
  --   tol       tolerance for relative error on two vectors;
  --   j         must be equal to sa'first.

  -- ON RETURN :
  --   j         the entry for which sa(j) equals x.

  function Number ( sols : Solution_List; flag : natural ) return natural;

  -- DESCRIPTION :
  --   Returns the number of solutions in the list with a multiplicity
  --   equal to flag.

  function Is_In ( sols : Solution_List; s : Solution; tol : double_float ) 
		 return boolean;
  function Is_In ( sa : Solution_Array; s : Solution; tol : double_float ) 
		 return boolean;

  -- DESCRIPTION :
  --   Returns true if the solution s belongs to the list or to the array.

  function Get ( sols : Solution_List; pos : positive ) return Solution;

  -- DESCRIPTION :
  --   Returns the solution at the given position.

  -- REQUIRED : pos <= Length_Of(sols).

-- CONSTRUCTORS :

  procedure Copy ( s1 : in Solution_List; s2 : in out Solution_List );
  procedure Copy ( s1 : in Solution_Array; s2 : in out Solution_Array );

  -- DESCRIPTION :
  --   Makes a deep copy of the list or the array of solutions.

  procedure Append ( first,last : in out Solution_List; s : in Solution );

  -- DESCRIPTION :
  --   The solution sol is appended to the list first;
  --   last is a pointer to the last element of the list first.

  procedure Add ( sols : in out Solution_List; s : in Solution );

  -- DESCRIPTION :
  --   The solution sol is appended to the list sols.

  procedure Add ( sols : in out Solution_List; s : in Solution;
                  tol : in double_float; other : out natural );

  -- DESCRIPTION :
  --   Append the solution to the list, if it does not already belong to it.

  procedure Change ( sols : in out Solution_List; pos : in positive;
                     s : in Solution; tol : in double_float; 
                     other : out natural );

  -- DESCRIPTION :
  --   Changes the solution at the given position into s, if the solution
  --   does not already occur.

  -- REQUIRED : pos <= Length_Of(sols).

  procedure Set_Continuation_Parameter
                ( sols : in out Solution_List; t : in double_complex );

  -- DESCRIPTION :
  --   All solutions in the list will be given the continuation parameter t.

  procedure Change_Multiplicity
                ( sols : in out Solution_List; pos : in positive;
                  m : in natural );

  -- DESCRIPTION :
  --   Changes the multiplicity of the solution with the given position
  --   into m.

  -- REQUIRED : pos <= Length_Of(sols).

  procedure Remove ( sols : in out Solution_List; pos : in positive );

  -- DESCRIPTION :
  --   Removes the solution with given position from the list.

  -- REQUIRED : pos <= Length_Of(sols).

  generic
    with function To_Be_Removed ( flag : in natural ) return boolean;
  procedure Delete ( sols : in out Solution_List );

  -- DESCRIPTION :
  --   Removes all solutions in s for which To_Be_Removed(s.m) holds.

  procedure Remove_All ( sols : in out Solution_List; flag : in natural );

  -- DESCRIPTION :
  --   All solutions with a multiplicity equal to flag are removed.

-- DESTRUCTORS :

  procedure Clear ( sa : in out Solution_Array );
  procedure Clear ( ls : in out Link_to_Solution );
  procedure Shallow_Clear ( sl : in out Solution_List );
  procedure    Deep_Clear ( sl : in out Solution_List );

  -- DESCRIPTION :
  --   Deallocation of the occupied memory space.
  --   A shallow clear only deallocates the pointers,
  --   so that the data may still be accessible by sharing,
  --   whereas a deep clear also makes the data inaccessible.

end Solutions;
