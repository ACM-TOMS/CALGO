with Integer_Vectors;   use Integer_Vectors;

package Integer_Graded_Lexicographical_Ordening is

-- DESCRIPTION :
--   This package provides a graded lexicographical ordening.

  function "<" ( v1,v2 : Vector ) return boolean;       -- return v1 < v2
  function "<" ( v1,v2 : Link_to_Vector ) return boolean;

  function ">" ( v1,v2 : Vector ) return boolean;       -- return v1 > v2
  function ">" ( v1,v2 : Link_to_Vector ) return boolean; 

end Integer_Graded_Lexicographical_Ordening;
