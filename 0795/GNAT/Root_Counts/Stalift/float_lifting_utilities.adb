package body Float_Lifting_Utilities is

  function Adaptive_Lifting ( l : Array_of_Lists ) return Vector is

    res : Vector(l'range);
    fac : constant double_float := 3.0;     -- multiplication factor
    max : constant double_float := 23.0;    -- upper bound for lifting

  begin
    for i in l'range loop
      res(i) := fac*double_float(Length_Of(l(i)));
      if res(i) > max
       then res(i) := max;
      end if;
    end loop;
    return res;
  end Adaptive_Lifting;

  procedure Search_Lifting ( l : in List; pt : in Vector;
                             found : out boolean; lif : out double_float ) is

    tmp : List := l;
    lpt : Link_to_Vector;

  begin
    found := false;
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt)
       then found := true;
            lif := lpt(lpt'last);
            exit;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
  end Search_Lifting;

  function Search_and_Lift ( l : List; pt : Vector ) return Vector is

    tmp : List := l;
    lpt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt)
       then return lpt.all;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return pt;
  end Search_and_Lift;

  function Search_and_Lift ( mic : Mixed_Cell; k : natural; pt : Vector )
                           return Vector is
  begin
    return Search_and_Lift(mic.pts(k),pt);
  end Search_and_Lift;

  function Induced_Lifting ( mixsub : Mixed_Subdivision; k : natural;
                             pt : Vector ) return Vector is

    tmp : Mixed_Subdivision := mixsub;
    res : Vector(pt'first..pt'last+1);

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        lpt : constant Vector := Search_and_Lift(mic,k,pt);
      begin
        if lpt'length = pt'length+1
         then return lpt;
         else tmp := Tail_Of(tmp);
        end if;
      end;
    end loop;
    res(pt'range) := pt;
    res(res'last) := 1.0;
    res(res'last) := Conservative_Lifting(mixsub,k,res);
    return res;
  end Induced_Lifting;

  function Induced_Lifting
               ( n : natural; mix : Integer_Vectors.Vector;
                 points : Array_of_Lists; mixsub : Mixed_Subdivision )
               return Array_of_Lists is

    res,res_last : Array_of_Lists(mix'range);
    cnt : natural := 1;
    tmp : List;

  begin
    for k in mix'range loop
      res_last(k) := res(k);
      tmp := points(cnt);
      while not Is_Null(tmp) loop
        declare
          pt : Link_to_Vector := Head_Of(tmp);
          lpt : constant Vector := Induced_Lifting(mixsub,k,pt.all);
        begin
          Append(res(k),res_last(k),lpt);
        end;
        tmp := Tail_Of(tmp);
      end loop;
      cnt := cnt + mix(k);
    end loop;
    return res;
  end Induced_Lifting;

  function Conservative_Lifting
               ( mic : Mixed_Cell; k : natural; point : Vector )
               return double_float is

    sp : double_float := mic.nor*Head_Of(mic.pts(k));
    spp : double_float:= mic.nor.all*point;
    res : double_float;

  begin
    if sp < spp
     then return point(point'last);
     else if mic.nor(mic.nor'last) = 0.0
           then res := point(point'last);
           else spp := spp - point(point'last)*mic.nor(mic.nor'last);
                res := (sp - spp)/mic.nor(mic.nor'last) + 1.0;
          end if;
          return res;
    end if;
  end Conservative_Lifting;

  function Conservative_Lifting ( mixsub : Mixed_Subdivision; k : natural;
                                  point : Vector ) return double_float is

    tmp : Mixed_Subdivision := mixsub;
    pt : Vector(point'range) := point;
    res : double_float;

  begin
    while not Is_Null(tmp) loop
      pt(pt'last) := Conservative_Lifting(Head_Of(tmp),k,pt);
      tmp := Tail_Of(tmp);
    end loop;
    res := pt(pt'last);
    Clear(pt);
    return res;
  end Conservative_Lifting;

end Float_Lifting_Utilities;
