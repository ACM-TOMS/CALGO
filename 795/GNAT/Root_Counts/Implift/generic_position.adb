with Lists_of_Integer_Vectors;         use Lists_of_Integer_Vectors;
with Power_Lists;                      use Power_Lists;
with Integer_Support_Functions;        use Integer_Support_Functions;

function Generic_Position ( p : Poly_Sys; tv : Tree_of_Vectors )
                          return boolean is

-- ALGORITHM :
--   Because tv contains the useful directions for the mixed volume,
--   only the cardinality of the faces of the first support list will be 
--   checked, for all top directions in tv.

  res : boolean := true;
  tmp : Tree_of_Vectors := tv;
  l : List := Construct_Power_List(p(p'first));

begin
  while not Is_Null(tmp) loop
    declare
      nd : Node := Head_Of(tmp);
      fc : List;
    begin
      exit when nd.d'length < Head_Of(l)'length;
      fc := Outer_Face(l,nd.d.all);
      res := (Length_Of(fc) <= 1);
      Deep_Clear(fc);
    end;
    exit when not res;
    tmp := Tail_Of(tmp);
  end loop;
  return res;
end Generic_Position;
