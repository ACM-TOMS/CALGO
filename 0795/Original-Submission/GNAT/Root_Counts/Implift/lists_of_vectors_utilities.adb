with Integer_Vectors_Utilities;       use Integer_Vectors_Utilities;
with Integer_Matrices;                use Integer_Matrices;
with Integer_Linear_System_Solvers;   use Integer_Linear_System_Solvers;  

package body Lists_of_Vectors_Utilities is

  procedure Compute_Normal ( v : in Integer_Vectors_of_Vectors.Vector;
                             n : out Link_to_Vector; deg : out natural ) is

    d : Link_to_Vector renames v(v'last);
    im : matrix(d'range,d'range);
    res : Link_to_Vector;
    cnt : integer;

  begin
    res := new Vector(d'range);
    cnt := im'first(1);
    for i in v'first..(v'last-1) loop
      for j in im'range(2) loop
        im(cnt,j) := v(i)(j) - d(j);
      end loop;
      cnt := cnt + 1;
    end loop;
    for i in cnt..im'last(1) loop
      for j in im'range(2) loop
        im(i,j) := 0;
      end loop;
    end loop;
    Upper_Triangulate(im);
    cnt := 1;
    for k in im'first(1)..im'last(1)-1 loop
      cnt := cnt*im(k,k);
    end loop;
    if cnt < 0
     then deg := -cnt;
     else deg := cnt;
    end if;
    Scale(im);
    Solve0(im,res.all);
    Normalize(res);
    n := res;
  end Compute_Normal;

  function Compute_Normal ( v : Integer_Vectors_of_Vectors.Vector )
                          return Link_to_Vector is

    deg : natural;
    res : Link_to_Vector;

  begin
    Compute_Normal(v,res,deg);
    return res;
  end Compute_Normal;

  function Pointer_to_Last ( l : List ) return List is

    res : List := l;

  begin
    if not Is_Null(res)
     then while not Is_Null(Tail_Of(res)) loop
            res := Tail_Of(res);
          end loop;
    end if;
    return res;
  end Pointer_to_Last;

  procedure Move_to_Front ( l : in out List; v : in Vector ) is

    tmp : List := l;
    found : boolean := false;
    first,lv : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lv := Head_Of(tmp);
      if Equal(lv.all,v)
       then found := true;
       else tmp := Tail_Of(tmp);
      end if;
      exit when found;
    end loop;
    if found
     then
       first := Head_Of(l);
       if first /= lv
        then
          lv.all := first.all;  Set_Head(tmp,lv);
          first.all := v;       Set_Head(l,first);
       end if;
    end if;
  end Move_to_Front;

  function Difference ( l1,l2 : List ) return List is

    res,res_last : List;
    tmp : List := l1;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if not Is_In(l2,pt.all)
       then Append(res,res_last,pt.all);
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Difference;

  function Different_Points ( l : List ) return List is

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      Append_Diff(res,res_last,Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Different_Points;

  procedure Remove_Duplicates ( l : in out List ) is

    res : List := Different_Points(l);

  begin
    Deep_Clear(l);
    l := res;
  end Remove_Duplicates;

end Lists_of_Vectors_Utilities;
