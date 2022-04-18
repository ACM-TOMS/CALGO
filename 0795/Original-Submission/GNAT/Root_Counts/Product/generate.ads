generic

  type boolean_array is array ( integer range <> ) of boolean;
  with procedure Process ( ar : in boolean_array; continue : out boolean );

procedure Generate ( k,first,last : in integer );

-- DESCRIPTION :
--   This procedure generates all possible unions of k
--   elements in the range first..last.
