with text_io,Natural_Vectors;       use text_io;
with Lists_of_Integer_Vectors;      use Lists_of_Integer_Vectors;

package Set_Structure is

-- DESCRIPTION :
--   This package manages the set structure of a polynomial system.

-- CONSTRUCTORS :

  procedure Init ( ns : in Natural_Vectors.Vector );

  -- DESCRIPTION :
  --   Initialization of the datastructures, must be the first operation.
  --   ns(i) contains the number of sets for the ith equation.

  procedure Add ( i,j,k : in natural );

  -- DESCRIPTION :
  --   To the jth set of the ith equation, the kth unknown will be added.

  procedure Remove ( i,j,k : in natural );

  -- DESCRIPTION :
  --   The kthe unknown will be removed from the jth set
  --   of the ith equation.

--  SELECTORS :

  function Empty return boolean;

  -- DESCRIPTION :
  --   Returns true if the set structure is empty.

  function Dimension return natural;

  -- DESCRIPTION :
  --   Returns the dimension, that is the number of equations.

  function Number_of_Sets ( i : natural ) return natural;

  -- DESCRIPTION :
  --   Returns the number of sets for the i-th equation.

  function Is_In ( i,j,k : natural ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the k-th unknown belongs to the j-th set of the
  --   i-th equation.

-- COMPUTING THE UPPER BOUND :

  function B return natural;

  -- REQUIRED :
  --   There exists a set structure for a polynomial system

  -- DESCRIPTION :
  --   Returns an upper bound to the number of solutions to
  --   a given polynomial system.

  procedure B ( bn : out natural; l : in out List );

  -- REQUIRED :
  --   There exists a set structure for a polynomial system

  -- DESCRIPTION :
  --   Returns the upper bound bn, together with the positions
  --   corresponding to the acceptable classes.

-- DESTRUCTOR :

  procedure Clear;

  -- DESCRIPTION :
  --   All used memory will be freed

end Set_Structure;
