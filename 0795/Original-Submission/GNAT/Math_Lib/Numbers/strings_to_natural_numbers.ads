package Strings_to_Natural_Numbers is

-- DESCRIPTION :
--   Natural numbers are extracted from strings, when possible.

  function Convert ( c : character ) return natural;

  -- DESCRIPTION :
  --   This function converts a character into a numerical value.
  --   If the character does not represent a decimal number, 
  --   then 10 is returned.

  function Convert ( s : string ) return integer;

  -- DESCRIPTION :
  --   This function scans the string, skips character that are
  --   no digits and converts the natural number in the string.
  --   If the string contains no natural number, then -1 is returned.

end Strings_to_Natural_Numbers;
