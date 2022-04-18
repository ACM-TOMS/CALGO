with Integer_Graded_Lexicographical_Ordening;
 use Integer_Graded_Lexicographical_Ordening;

package body Integer_Support_Functions is

  function Maximal_Support ( l : List; v : Vector ) return integer is

    sp,max : integer;
    tmp : List;

  begin
    if not Is_Null(l)
     then max := Head_Of(l).all*v;
          tmp := Tail_Of(l);
          while not Is_Null(tmp) loop
            sp := Head_Of(tmp).all*v;
            if sp > max
             then max := sp;
            end if;
            tmp := Tail_Of(tmp);
          end loop;
          return max;
     else return 0;
    end if;
  end Maximal_Support;

  function Minimal_Support ( l : List; v : Vector ) return integer is

    sp,min : integer;
    tmp : List;

  begin
    if not Is_Null(l)
     then min := Head_Of(l).all*v;
          tmp := Tail_Of(l);
          while not Is_Null(tmp) loop
            sp := Head_Of(tmp).all*v;
            if sp < min
             then min := sp;
            end if;
            tmp := Tail_Of(tmp);
          end loop;
          return min;
     else return 0;
    end if;
  end Minimal_Support;

  procedure Min_Max ( l : in List; k : in integer;
                      min,max : in out integer ) is

    tmp : List;
    v : Link_to_Vector;

  begin
    if not Is_Null(l)
     then tmp := l;
          v := Head_Of(tmp);
          min := v(k);  max := min;
          tmp := Tail_Of(tmp);
          while not Is_Null(tmp) loop
            v := Head_Of(tmp);
            if v(k) < min
             then min := v(k);
             elsif v(k) > max
                 then max := v(k);
            end if;
            tmp := Tail_Of(tmp);
          end loop;
    end if;
  end Min_Max;

  function Graded_Max ( l : List ) return Link_to_Vector is

    res : Link_to_Vector := new Vector'(Head_Of(l).all);
    tmp : List := Tail_Of(l);
    ele : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      ele := Head_Of(tmp);
      if ele > res
       then res.all := ele.all;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Graded_Max;

  function Face ( l : List; v : Vector; m : integer ) return List is

    res,tmp,res_last : List;
    d : Vector(v'range);

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      d := Head_Of(tmp).all;
      if d*v = m
       then Append(res,res_last,d);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Face;

  function Inner_Face ( l : List; v : Vector ) return List is
  begin
    return Face(l,v,Minimal_Support(l,v));
  end Inner_Face;

  function Outer_Face ( l : List; v : Vector ) return List is
  begin
    return Face(l,v,Maximal_Support(l,v));
  end Outer_Face;

end Integer_Support_Functions;
