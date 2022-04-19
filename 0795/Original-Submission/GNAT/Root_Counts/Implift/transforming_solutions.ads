with Transformations,Solutions;          use Transformations,Solutions;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers,Complex_Vectors;

package Transforming_Solutions is

-- DESCRIPTION :
--   This package offers some routines for transforming solutions.

  procedure Transform ( t : in Transfo; s : in out Solution );
  function  Transform ( t : Transfo; s : Solution ) return Solution;

  procedure Transform ( t : in Transfo; l : in out Solution_List );
  function  Transform ( t : Transfo; l : Solution_List ) return Solution_List;

  -- DESCRIPTION :
  --   The transformation t will be applied.

  function Insert ( c : double_complex; i : integer; s : Solution )
                  return Solution;

  procedure Insert ( c : in double_complex; i : in integer; 
		     l : in out Solution_List );
  function Insert ( c : double_complex; i : integer; l : Solution_List )
		  return Solution_List;

  -- DESCRIPTION :
  --   To the ith component of the solution vector, c will be inserted.

  function Insert ( cv : vector; i : integer; s : Solution )
		  return Solution_List;

  -- DESCRIPTION :
  --   All components in the vector cv will be inserted as the
  --   the ith component in the solution vector.

end Transforming_Solutions;
