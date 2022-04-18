package body Transforming_Solutions is

  procedure Transform ( t : in Transfo; s : in out Solution ) is
  begin
    Apply(t,s.v);
  end Transform;

  function  Transform ( t : Transfo; s : Solution ) return Solution is

    res : Solution(s.n);

  begin
    res.m := s.m;
    res.t := s.t;
    res.v := t*s.v;
    return res;
  end Transform;

  procedure Transform ( t : in Transfo; l : in out Solution_List ) is

    tmp : Solution_List;

  begin
    if not Is_Null(l)
     then declare
	    n : natural := Head_Of(l).n;
	    s : Solution(n);
          begin
	    tmp := l;
            while not Is_Null(tmp) loop
              Apply(t,Head_Of(tmp).v);
	      tmp := Tail_Of(tmp);
            end loop;
          end;
    end if;
  end Transform;

  function Transform ( t : Transfo; l : Solution_List ) return Solution_List is

    res,res_last,tmp : Solution_List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      Append(res,res_last,Transform(t,Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Transform;

  function Insert ( c : double_complex; i : integer; s : Solution )
                  return Solution is

    res : Solution(s.n+1);

  begin
    res.m := s.m;
    res.t := s.t;
    for j in res.v'first..(i-1) loop
      res.v(j) := s.v(j);
    end loop;
    res.v(i) := c;
    for j in (i+1)..res.v'last loop
      res.v(j) := s.v(j-1);
    end loop;
    return res;
  end Insert;

  procedure Insert ( c : in double_complex; i : in integer;
		     l : in out Solution_List ) is
  begin
    if not Is_Null(l)
     then declare
	    tmp : Solution_List;
	    res,res_last : Solution_List;
	    sol : Solution(Head_Of(l).n+1);
          begin
	    tmp := l;
	    while not Is_Null(tmp) loop
	      declare
		ls : Link_to_Solution := Head_Of(tmp);
              begin
		sol.m := ls.m;
		sol.t := ls.t;
		for j in sol.v'first..(i-1) loop
		  sol.v(j) := ls.v(j);
                end loop;
		sol.v(i) := c;
		for j in (i+1)..sol.v'last loop
		  sol.v(j) := ls.v(j-1);
                end loop;
		Append(res,res_last,sol);
              end;
	      tmp := Tail_Of(tmp);
            end loop;
	    Clear(l); l := res;
          end;
    end if;
  end Insert;

  function Insert ( c : double_complex; i : integer; l : Solution_List )
		  return Solution_List is

    res : Solution_List;

  begin
    if not Is_Null(l)
     then declare 
	    tmp,res_last : Solution_List;
          begin
	    tmp := l;
	    while not Is_Null(tmp) loop
	      Append(res,res_last,Insert(c,i,Head_Of(tmp).all));
              tmp := Tail_Of(tmp);
            end loop;
	  end;
    end if;
    return res;
  end Insert;

  function Insert ( cv : Vector; i : integer; s : Solution )
		  return Solution_List is

    res,res_last : Solution_List;

  begin
    for j in cv'range loop
      Append(res,res_last,Insert(cv(j),i,s));
    end loop;
    return res; 
  end Insert;

end Transforming_Solutions;
