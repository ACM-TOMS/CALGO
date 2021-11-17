with Integer_Vectors;               use Integer_Vectors;
with Integer_Vectors_of_Vectors;

package body Common_Faces_of_Polytope is

  function Have_Common_Point
              ( l : List; v : Integer_Vectors_of_Vectors.Vector ) 
              return boolean is

  -- DESCRIPTION :
  --   Returns true if at least one point in v belongs to the list l.

  begin
    for i in v'range loop
      if Is_In(l,v(i).all)
       then return true;
      end if;
    end loop;
    return false;
  end Have_Common_Point;

  function Is_Neighbor1 ( l : List; fc : Face ) return boolean is
  begin
    return Have_Common_Point(l,fc.all);
  end Is_Neighbor1;

  function Is_Neighbor ( l : List; fc : Face ) return boolean is

    cntnotin : natural := 0;
       -- counts the points in the face fc that are not in the list l

  begin
    for i in fc'range loop
      if not Is_In(l,fc(i).all)
       then cntnotin := cntnotin + 1;
            if cntnotin > 1
             then return false;
            end if;
      end if;
    end loop;
    return true;
  end Is_Neighbor;

  function Neighboring_Faces ( mic : Mixed_Cell; fs : Faces; i : natural )
                             return Faces is

    tmp : Faces := fs;
    res,res_last : Faces;

  begin
    while not Is_Null(tmp) loop
      declare
        fc : Face := Head_Of(tmp);
      begin
        if Is_Neighbor(mic.pts(i),fc)
         then Append(res,res_last,fc);
        end if;
        tmp := Tail_Of(tmp);
      end;
    end loop;
    return res;
  end Neighboring_Faces;

  function Neighboring_Faces ( mic : Mixed_Cell; afs : Array_of_Faces )
                             return Array_of_Faces is

    res : Array_of_Faces(afs'range);

  begin
    for i in res'range loop
      res(i) := Neighboring_Faces(mic,afs(i),i);
    end loop;
    return res;
  end Neighboring_Faces;

end Common_Faces_of_Polytope;
