with Natural_Vectors;   use Natural_Vectors;

package Natural_Graded_Lexicographical_Ordening is

-- DESCRIPTION :
--   This package provides a graded lexicographical ordening.

  function "<" ( v1,v2 : Vector ) return boolean;       -- return v1 < v2
  function "<" ( v1,v2 : Link_to_Vector ) return boolean;

  function ">" ( v1,v2 : Vector ) return boolean;       -- return v1 > v2
  function ">" ( v1,v2 : Link_to_Vector ) return boolean; 

end Natural_Graded_Lexicographical_Ordening;
