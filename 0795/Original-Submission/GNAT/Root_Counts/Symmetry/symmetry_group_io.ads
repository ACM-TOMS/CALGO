with text_io;                      use text_io;
with Permutations,Symmetry_Group;  use Permutations,Symmetry_Group;

package Symmetry_Group_io is

-- DESCRIPTION :
--   This package contains routines for input and output of 
--   permutations and groups of permutations.

  procedure get ( p : out Permutation );
  procedure get ( file : in file_type; p : out Permutation );

  -- DESCRIPTION :
  --   Reads a permutation from standard input or from file.
   
  procedure get ( l : in out List_of_Permutations; n,nb : in natural );
  procedure get ( file : in file_type;
		  l : in out List_of_Permutations; n,nb : in natural );

  -- DESCRIPTION :
  --   Reads a list of permutations from standard in put or from file.
  --   Vectors that do not represent permutations are ignored.

  -- ON ENTRY :
  --   n          the number of elements in the permutations;
  --   nb         the total number of permutations that must be read.
      
  procedure put ( p : in Permutation );
  procedure put ( file : in file_type; p : in Permutation );

  -- DESCRIPTION :
  --   Writes a permutation on standard output or on file.

  procedure put ( l : in List_of_Permutations );
  procedure put ( file : in file_type; l : in List_of_Permutations );

  -- DESCRIPTION :
  --   Writes a list of permutations on standard output or on file.

end Symmetry_Group_io;
