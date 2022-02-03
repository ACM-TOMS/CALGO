package body Lists_of_Float_Vectors is

-- CONSTRUCTORS :

  function Deep_Create ( v : Float_Vectors_of_Vectors.Vector ) return List is

    res,res_last : List;

  begin
    for i in v'range loop
      Append(res,res_last,v(i).all);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( v : Float_Vectors_of_Vectors.Vector )
                          return List is

    res,res_last : List;

  begin
    for i in v'range loop
      Append(res,res_last,v(i));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( l : List ) return Float_Vectors_of_Vectors.Vector is

    res : Float_Vectors_of_Vectors.Vector(1..Length_Of(l));
    tmp : List := l;

  begin
    for i in res'range loop
      res(i) := new vector'(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( l : List )
                          return Float_Vectors_of_Vectors.Vector is

    res : Float_Vectors_of_Vectors.Vector(1..Length_Of(l));
    tmp : List := l;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shallow_Create;

  procedure Copy ( l1 : in List; l2 : in out List ) is
    tmp,l2_last : List;
    lv : Link_to_Vector;
  begin
    Deep_Clear(l2);
    tmp := l1;
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      Append(l2,l2_last,lv.all);
      tmp := Tail_Of(tmp);
    end loop;
  end Copy;

  procedure Append ( first,last : in out List; v : in Vector ) is

    lv : Link_to_Vector := new Vector'(v);

  begin
    if Is_Null(first)
     then Construct(lv,first);
          last := first;
     else declare
            tmp : List;
          begin
            Construct(lv,tmp);
            Swap_Tail(last,tmp);
            last := Tail_Of(last);
          end;
    end if;
  end Append;

  procedure Append_Diff ( first,last : in out List; v : in Vector ) is
  begin
    if not Is_In(first,v)
     then Append(first,last,v);
    end if;
  end Append_Diff;

  procedure Append_Diff ( first,last : in out List; v : in Link_to_Vector ) is
  begin
    if v /= null and then not Is_In(first,v)
     then Append(first,last,v);
    end if;
  end Append_Diff;

  procedure Deep_Concat ( first,last : in out List; l : in List ) is

    tmp : List;
    lv : Link_to_Vector;

  begin
    if not Is_Null(l)
     then tmp := l;
          while not Is_Null(tmp) loop
            lv := Head_Of(tmp);
            Append(first,last,lv.all);
            tmp := Tail_Of(tmp);
          end loop;
    end if;
  end Deep_Concat;

  procedure Shallow_Concat ( first,last : in out List; l : in List ) is
  begin
    Concat(first,last,l);
  end Shallow_Concat;

  procedure Deep_Concat_Diff ( first,last : in out List; l : in List ) is

    tmp : List;
    lv : Link_to_Vector;

  begin
    if not Is_Null(l)
     then tmp := l;
          while not Is_Null(tmp) loop
            lv := Head_Of(tmp);
            Append_Diff(first,last,lv.all);
            tmp := Tail_Of(tmp);
          end loop;
    end if;
  end Deep_Concat_Diff;

  procedure Shallow_Concat_Diff ( first,last : in out List; l : in List ) is

    tmp : List;
    lv : Link_to_Vector;

  begin
    if not Is_Null(l)
     then tmp := l;
          while not Is_Null(tmp) loop
            lv := Head_Of(tmp);
            Append_Diff(first,last,lv);
            tmp := Tail_Of(tmp);
          end loop;
    end if;
  end Shallow_Concat_Diff;

  procedure Remove ( l : in out List; x : in Vector ) is

    lpt : Link_to_Vector;
    found : boolean;
    l1,l2 : List;

  begin
    if not Is_Null(l)
     then
       lpt := Head_Of(l);
       if lpt.all = x
        then Clear(lpt);
             l := Tail_Of(l);
        else found := false;
             l1 := l;
             l2 := Tail_Of(l1);
             while not Is_Null(l2) loop
               lpt := Head_Of(l2);
               found := (lpt.all = x);
               exit when found;
               l1 := l2;
               l2 := Tail_Of(l1);
             end loop;
             if found
              then Clear(lpt);
                   l2 := Tail_Of(l2);
                   Swap_Tail(l1,l2);
             end if;
       end if;
    end if;
  end Remove;

  procedure Remove ( l : in out List; x : in Link_to_Vector ) is
  begin
    if x /= null
     then Remove(l,x.all);
    end if;
  end Remove;

  procedure Swap_to_Front ( l : in out List; x : in Vector ) is

    first : Link_to_Vector;
    pt : Link_to_Vector;
    tmp : List;
    done : boolean := false;

  begin
    if not Is_Null(l)
     then first := Head_Of(l);
          if first.all /= x
           then tmp := Tail_Of(l);
                while not Is_Null(tmp) loop
                  pt := Head_Of(tmp);
                  if pt.all = x
                   then Set_Head(tmp,first);
                        Set_Head(l,pt);
                        done := true;
                  end if;
                  exit when done;
                  tmp := Tail_Of(tmp);
                end loop;
          end if;
    end if;
  end Swap_to_Front;

  procedure Swap_to_Front ( l : in out List; x : in Link_to_Vector ) is
  begin
    if x /= null
     then Swap_to_Front(l,x.all);
    end if;
  end Swap_to_Front;

-- SELECTORS :

  function Is_In ( l : List; v : Vector ) return boolean is

    tmp : List;
    v2 : Link_to_Vector;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      v2 := Head_Of(tmp);
      if Equal(v2.all,v)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( l : List; v : Link_to_Vector ) return boolean is
  begin
    if v = null
     then return false;
     else return Is_In(l,v.all);
    end if;
  end Is_In;

  function Is_Equal ( l1,l2 : List ) return boolean is

    function Check_All ( l1,l2 : List ) return boolean is

     -- DESCRIPTION :
     --   Returns true if all elements of l1 occur in l2.

      tmp : List;
      v : Link_to_Vector;

    begin
      tmp := l1;
      while not Is_Null(tmp) loop
	v := Head_Of(tmp);
	if not Is_In(l2,v)
	 then return false;
	 else tmp := Tail_Of(tmp);
        end if;
      end loop;
      return true;
    end Check_All;

  begin
    if not Check_All(l1,l2)
     then return false;
     elsif not Check_All(l2,l1)
	 then return false;
	 else return true;
    end if;
  end Is_Equal;

-- DESTRUCTORS :

  procedure Deep_Clear ( l : in out List ) is

    tmp : List;
    v : Link_to_Vector;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      v := Head_Of(tmp);
      Clear(v);
      tmp := Tail_Of(tmp);
    end loop;
    Shallow_Clear(l);
  end Deep_Clear;

  procedure Shallow_Clear ( l : in out List ) is
  begin
    Lists_of_Link_to_Float_Vectors.Clear
                       (Lists_of_Link_to_Float_Vectors.List(l));
  end Shallow_Clear;

end Lists_of_Float_Vectors;
