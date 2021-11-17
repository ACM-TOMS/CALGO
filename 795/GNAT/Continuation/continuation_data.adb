package body Continuation_Data is

-- AUXILIARY :

  procedure Append ( first,last : in out Solution_List;
                     ls : in Link_to_Solution ) is
  begin
    if Is_Null(first)
     then Construct(ls,first);
          last := first;
     else declare
            tmp : Solution_List;
          begin
            Construct(ls,tmp);
            Swap_Tail(last,tmp);
            last := Tail_Of(last);
          end;
    end if;
  end Append;

-- CREATORS :

  function Shallow_Create ( s : Link_to_Solution ) return Solu_Info is

    res : Solu_Info;

  begin
    res.sol := s;
    Init_Info(res);
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solution ) return Solu_Info is

    res : Solu_Info;

  begin
    res.sol := new Solution'(s);
    Init_Info(res);
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solution_Array ) return Solu_Info_Array is

    res : Solu_Info_Array(s'range);

  begin
    for k in res'range loop
      res(k) := Shallow_Create(s(k));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solution_Array ) return Solu_Info_Array is

    res : Solu_Info_Array(s'range);

  begin
    for k in res'range loop
      res(k) := Deep_Create(s(k).all);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solution_List )  return Solu_Info_Array is
  
    res : Solu_Info_Array(1..Length_Of(s));
    tmp : Solution_List := s;

  begin
    for k in res'range loop
      res(k) := Shallow_Create(Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shallow_Create;
   
  function Deep_Create ( s : Solution_List )  return Solu_Info_Array is
   
    res : Solu_Info_Array(1..Length_Of(s));
    tmp : Solution_List := s;
 
  begin
    for k in res'range loop
      res(k) := Deep_Create(Head_Of(tmp).all);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solu_Info ) return Link_to_Solution is
  begin
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
    return s.sol;
  end Shallow_Create;

  function Deep_Create ( s : Solu_Info ) return Solution is
  begin
    s.sol.err := s.cora;
    s.sol.rco := s.rcond;
    s.sol.res := s.resa;
    return s.sol.all;
  end Deep_Create;

  function Shallow_Create ( s : Solu_Info_Array ) return Solution_Array is

    res : Solution_Array(s'range);

  begin
    for k in s'range loop
      res(k) := Shallow_Create(s(k));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solu_Info_Array ) return Solution_Array is

    res : Solution_Array(s'range);

  begin
    for k in s'range loop
      res(k) := new Solution'(Deep_Create(s(k)));
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( s : Solu_Info_Array ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for k in s'range loop
      Append(res,res_last,Shallow_Create(s(k)));
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( s : Solu_Info_Array ) return Solution_List is

    res,res_last : Solution_List;

  begin
    for k in s'range loop
      Append(res,res_last,Deep_Create(s(k)));
    end loop;
    return res;
  end Deep_Create;

-- OPERATIONS ON Solu_Info :

  procedure Copy_Info ( s1 : in Solu_Info; s2 : in out Solu_Info ) is
  begin
    s2.corr := s1.corr; s2.cora := s1.cora;
    s2.resr := s1.resr; s2.resa := s1.resa;
    s2.rcond := s1.rcond;  s2.length_path := s1.length_path;
    s2.nstep := s1.nstep; s2.nfail := s1.nfail;
    s2.niter := s1.niter; s2.nsyst := s1.nsyst;
  end Copy_Info;

  procedure Copy_Solu ( s1 : in Solu_Info; s2 : in out Solu_Info ) is
  begin
    Clear(s2.sol);
    s2.sol := new Solution'(s1.sol.all);
  end Copy_Solu;

  procedure Copy ( s1 : in Solu_Info; s2 : in out Solu_Info ) is
  begin
    Copy_Info(s1,s2);
    Copy_Solu(s1,s2);
  end Copy;

  procedure Init_Info ( s : in out Solu_Info ) is
  begin
    s.corr := 0.0; s.cora := 0.0; s.resr := 0.0; s.resa := 0.0; s.rcond := 0.0;
    s.length_path := 0.0;
    s.nstep := 0; s.nfail := 0; s.niter := 0; s.nsyst := 0;
  end Init_Info;

  procedure Add_Info ( s1 : in out Solu_Info; s2 : in Solu_Info ) is
  begin
    s1.nstep := s1.nstep + s2.nstep;
    s1.nfail := s1.nfail + s2.nfail;
    s1.niter := s1.niter + s2.niter;
    s1.nsyst := s1.nsyst + s2.niter;
    s1.length_path := s1.length_path + s2.length_path;
  end Add_Info;

  procedure Update_Info ( s1 : in out Solu_Info; s2 : in Solu_Info ) is
  begin
    s1.corr := s2.corr; s1.cora := s2.cora;
    s1.resr := s2.resr; s1.resa := s2.resa;
    s1.rcond := s2.rcond;
    Add_Info(s1,s2);
  end Update_Info;

-- OPERATIONS ON Solu_Info_Array :

  procedure Copy ( s : in Solu_Info_Array; sa : in out Solution_Array ) is
  begin
    Clear(sa);
    for k in s'range loop
      sa(k) := new Solution'(s(k).sol.all);
    end loop;
  end Copy;

  procedure Copy ( sa : in Solution_Array; s : in out Solu_Info_Array ) is
  begin
    for k in sa'range loop
      Clear(s(k).sol);
      s(k).sol := new Solution'(sa(k).all);
    end loop;
  end Copy;

-- DESTRUCTORS :

  procedure Clear ( s : in out Solu_Info ) is
  begin
    Clear(s.sol);
  end Clear;

  procedure Clear ( s : in out Solu_Info_Array ) is
  begin
    for k in s'range loop
      Clear(s);
    end loop;
  end Clear;

end Continuation_Data;
