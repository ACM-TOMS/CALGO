with Sets_of_Unknowns;          use Sets_of_Unknowns;

package Partitions_of_Sets_of_Unknowns is

-- DESCRIPTION :
--   This package provides a data abstraction for enumerating all
--   partitions of a given set of unknowns.

  type partition is array ( positive range <> ) of Set;

-- CREATORS :

  procedure Create ( p : in out Partition; n : in natural );

  -- DESCRIPTION :
  --   Creates all sets in the partition, to be ready to contain
  --   at most n unknowns.

  function Create ( p : Partition ) return Partition;
 
  -- DESCRIPTION :
  --   Returns a new partition which is an exact copy of the given one.

-- CONSTRUCTORS :

  generic 
    with procedure Process ( p : in Partition; continue : out boolean );
  procedure Generate_Partitions ( s : in Set );

  -- DESCRIPTION :
  --   Generates all partitions of a given set of unknowns.
  --   The procedure Process is invoked each time a new partition is
  --   generated.  The generation can be stopped by setting
  --   continue to false.  Otherwise, when continue is set to true,
  --   the generation continues.

  -- NOTE :
  --   While processing the partition, it might be needed to copy
  --   the resulting partition, as sharing occurs.

-- SELECTOR :

  function Number_of_Partitions ( n : natural ) return natural;

  -- DESCRIPTION :
  --   Returns the number of partitions of a set of n unknowns.

-- DESTRUCTOR :

  procedure Clear ( p : in out Partition );

  -- DESCRIPTION :
  --   Deallocates the occupied memory.

end Partitions_of_Sets_of_Unknowns;
