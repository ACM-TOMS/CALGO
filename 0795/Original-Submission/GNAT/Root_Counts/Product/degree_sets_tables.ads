with Integer_Matrices;     use Integer_Matrices;
with Sets_of_Unknowns;     use Sets_of_Unknowns;

package Degree_Sets_Tables is

-- DESCRIPTION :
--   This package provides facilities for computing generalized permanents,
--   based on a set structure.

-- DATASTRUCTURES :

  type Array_of_Sets is array ( integer range <> ) of Set;
  
  type Degree_Sets_Table ( n,m : natural ) is record
    s : Array_of_Sets(1..m);
    a : Matrix(1..n,1..m);
  end record;

-- CONSTRUCTOR :

  function Create return Degree_Sets_Table;

  -- DESCRIPTION :
  --   Selects the information from the package Set_Structure to create
  --   a degree set structure table.

-- PERMANENT COMPUTATIONS :

  function Permanent ( dst : Degree_Sets_Table ) return natural;

  -- DESCRIPTION :
  --   Returns the generalized permanent, based on the set structure.

-- DESTRUCTOR :

  procedure Clear ( ase : in out Array_of_Sets );
  procedure Clear ( dst : in out Degree_Sets_Table );

  -- DESCRIPTION :
  --   Deallocates the space occupied by the array of sets.

end Degree_Sets_Tables;
