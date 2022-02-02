with Integer_Vectors_Utilities;       use Integer_Vectors_Utilities;

package body Transforming_Integer_Vector_Lists is

  procedure Shift ( l : in out List; v : in Vector ) is

    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      declare
        lv : Link_to_Vector := Head_Of(tmp);
      begin
        lv.all := lv.all - v;
        Set_Head(tmp,lv);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Shift;

  procedure Shift ( l : in out List; v : in Link_to_Vector ) is
  begin
    if v /= null
     then Shift(l,v.all);
    end if;
  end Shift;

  function Shift ( l : List; v : Vector ) return List is

    tmp,res,res_last : List;
    v1 : Vector(v'range);

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      v1 := Head_Of(tmp).all;
      declare
        v2 : Vector(v1'range) := v1 - v;
      begin
        Append(res,res_last,v2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shift;

  function Shift ( l : List; v : Link_to_Vector ) return List is
  begin
    if v = null
     then declare
            res : List;
          begin
            Copy(l,res);
            return res;
          end;
     else return Shift(l,v.all);
    end if;
  end Shift;

  function "*"( l : List; t : Transfo ) return List is
  begin
    return t*l;
  end "*";

  function "*"( t : Transfo; l : List ) return List is

    tmp,res,res_last : List;
    d,td : Link_to_Vector;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      d := Head_Of(tmp);
      td := t*d;
      Append(res,res_last,td);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end "*";

  procedure Apply ( l : in out List; t : in Transfo ) is

    res : List := t*l;

  begin
    Copy(res,l);
  end Apply;

  function Reduce ( l : List; i : integer ) return List is

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        d1 : Link_to_Vector := Head_Of(tmp);
        d2 : Link_to_Vector := Reduce(d1,i);
      begin
       -- Append_Diff(res,res_last,d2);      -- be aware of duplicates !
        Append(res,res_last,d2);      -- be aware of duplicates !
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Reduce;

  procedure Reduce ( l : in out List; i : in integer ) is

    res : List := Reduce(l,i);

  begin
    Copy(res,l);
  end Reduce;

  function Insert ( l : List; i,a : integer ) return List is

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        d1 : Link_to_Vector := Head_Of(tmp);
        d2 : Link_to_Vector := Insert(d1,i,a);
      begin
        Append(res,res_last,d2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert;

  procedure Insert ( l : in out List; i,a : in integer ) is

    res : List := Insert(l,i,a);

  begin
    Deep_Clear(l);
    l := res;
  end Insert;

  function Transform_and_Reduce ( t : Transfo; i : integer; l : List )
                                return List is
    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        d  : Link_to_Vector := Head_Of(tmp);
        td : Vector(d'range) := t*d.all;
        dr : Link_to_Vector := new Vector'(Reduce(td,i));
      begin
        Append(res,res_last,dr);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Transform_and_Reduce;

  procedure Transform_and_Reduce ( t : in Transfo; i : in integer;
                                   l : in out List ) is

    res : List := Transform_and_Reduce(t,i,l);

  begin
    Deep_Clear(l);
    l := res;
  end Transform_and_Reduce;

  function Insert_and_Transform ( l : List; i,a : integer; t : Transfo )
                                return List is

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        d : Link_to_Vector := Insert_and_Transform(Head_Of(tmp),i,a,t);
      begin
        Append(res,res_last,d);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Insert_and_Transform;

  procedure Insert_and_Transform
             ( l : in out List; i,a : in integer; t : in Transfo ) is

    res : List := Insert_and_Transform(l,i,a,t);

  begin
    Deep_Clear(l);
    l := res;
  end Insert_and_Transform;

end Transforming_Integer_Vector_Lists;
