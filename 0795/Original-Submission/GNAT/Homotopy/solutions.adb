with unchecked_deallocation;
with Floating_Equalities;        use Floating_Equalities;

package body Solutions is

  use List_of_Solutions;

-- CREATORS :

  function Create ( sl : Solution_List ) return Solution_Array is

    sa : Solution_Array(1..Length_Of(sl));

  begin
    if not Is_Null(sl)
     then declare
            i : positive := 1;
            temp : Solution_List := sl;
          begin
            while not Is_Null(temp) loop
              sa(i) := new Solution'(Head_Of(temp).all);
              i := i + 1;
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
    return sa;
  end Create;

  function Create ( sa : Solution_Array ) return Solution_List is

    sl : Solution_List;

  begin
    if sa'first <= sa'last
     then declare
            n : natural := sa(sa'first).n;
            sol : Solution(n) := sa(sa'first).all;
            l : Link_to_Solution := new Solution'(sol);
            last,tmp : Solution_List;
          begin
            Construct(l,sl);
            last := sl;
            for i in (sa'first+1)..sa'last loop
              sol := sa(i).all;
              l := new Solution'(sol);
              Construct(l,tmp);
              Swap_Tail(last,tmp);
              last := Tail_Of(last);
            end loop;
          end;
    end if;
    return sl;
  end Create;

-- SELECTORS :

  function Equal ( s1,s2 : Solution; tol : double_float ) return boolean is
  begin
    if (s1.t /= s2.t) or else (s1.n /= s2.n)
     then return false;
     else return Is_Equal(s1.v,s2.v,tol);
    end if;
  end Equal;

  function Equal ( s1,s2 : Solution_List; tol : double_float )
                 return boolean is
  begin
    if Is_Null(s1) and Is_Null(s2)
     then return true;
     elsif Is_Null(s1) or Is_Null(s2)
         then return false;
         else declare
                temp1 : Solution_List := s1;
                temp2 : Solution_List := s2;
              begin
                While not Is_Null(temp1) and not Is_Null(s2) loop
                  if not Equal(Head_Of(temp1).all,Head_Of(temp2).all,tol)
                   then return false;
                   else temp1 := Tail_Of(temp1);
                        temp2 := Tail_Of(temp2);
                  end if;
                end loop;
                if Is_Null(temp1) and Is_Null(temp2)
                 then return true;
                 else return false;
                end if;
              end;
    end if;
  end Equal;

  function Equal ( s1,s2 : Solution_Array; tol : double_float )
                 return boolean is
  begin
    if s1'first /= s2'first
     then return false;
     elsif s1'last /= s2'last
         then return false;
         else for i in s1'range loop
                if not Equal(s1(i).all,s2(i).all,tol)
                 then return false;
                end if;
              end loop;
    end if;
    return true;
  end Equal;

  procedure Equals ( sols : in out Solution_List; flag : in natural;
                     tol : in double_float; same : out boolean ) is
  begin
    same := false;
    if not Is_Null(sols)
     then declare
            n : natural := Head_Of(sols).n;
            i : natural := 1;
            s1,s2 : Solution(n);
            temp : Solution_List := sols;
          begin
            while not Is_Null(temp) loop
              s1 := Head_Of(temp).all;
              for j in (i+1)..Length_Of(sols) loop
                s2 := Get(sols,j);
                if Equal(s1,s2,tol)
                 then same := true;
                      Change_Multiplicity(sols,i,flag);
                      Change_Multiplicity(sols,j,flag);
                end if;
              end loop;
              temp := Tail_Of(temp);
              i := i + 1;
            end loop;
          end;
    end if;
  end Equals;

  procedure Equals ( sa : in Solution_Array; x : in Vector; i : in natural;
                     tol : in double_float; j : in out natural ) is

    eq : boolean;

  begin
    while j < i loop
      eq := true;
      for k in x'range loop
        if modulus(sa(j).v(k) - x(k)) > tol
         then eq := false;
        end if;
        exit when not eq;
      end loop;
      exit when eq;
      j := j + 1;
    end loop;
  end Equals;

  function Number ( sols : Solution_List; flag : natural ) return natural is

    res : natural := 0;

  begin
    if Is_Null(sols)
     then return res;
     else declare
            temp : Solution_List := sols;
            ls : Link_to_Solution;
          begin
            while not Is_Null(temp) loop
              if Head_Of(temp).m = flag
               then res := res + 1;
              end if;
              temp := Tail_Of(temp);
            end loop;
          end;
          return res;
    end if;
  end Number;

  function Is_In ( sols : Solution_List; s : Solution; tol : double_float )
                 return boolean is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      if Equal(Head_Of(tmp).all,s,tol)
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Is_In ( sa : Solution_Array; s : Solution; tol : double_float )
                 return boolean is
  begin
    for i in sa'range loop
      if Equal(sa(i).all,s,tol)
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

  function Get ( sols : Solution_List; pos : positive )
               return Solution is
  begin
    if pos <= Length_Of(sols)
     then declare
            temp : Solution_List := sols;
            count : natural := 1;
          begin
            while not Is_Null(temp) loop
              if count = pos
               then return Head_Of(temp).all;
               else temp := Tail_Of(temp);
                    count := count + 1;
              end if;
            end loop;
          end;
    end if;
    declare
      s : Solution(0);
    begin
      return s;
    end;
  end Get;

-- CONSTRUCTORS :

  procedure Copy ( s1 : in Solution_List; s2 : in out Solution_List ) is
  begin
    Clear(s2);
    if not Is_Null(s1)
     then declare
            temp : Solution_List := s1;
            last : Solution_List;
            n : natural := Head_Of(s1).n;
            sol : Solution(n) := Head_Of(temp).all;
          begin
            declare
              l : Link_to_Solution := new Solution'(sol);
            begin
              Construct(l,s2);
            end;
            last := s2;
            temp := Tail_Of(temp);
            while not Is_Null(temp) loop
              sol := Head_Of(temp).all;
              declare
                l : Link_to_Solution := new Solution'(sol);
                tmp : Solution_List;
              begin
                Construct(l,tmp);
                Swap_Tail(last,tmp);
              end;
              last := Tail_Of(last);
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
  end Copy;

  procedure Copy ( s1 : in Solution_Array; s2 : in out Solution_Array ) is
  begin
    Clear(s2);
    for i in s1'range loop
      s2(i) := new Solution'(s1(i).all);
    end loop;
  end Copy;

  procedure Append ( first,last : in out Solution_List; s : in Solution ) is

    ls : Link_to_Solution := new Solution'(s);

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

  procedure Add ( sols : in out Solution_List; s : in Solution ) is

    last,temp,tmp : Solution_List;
    ls : Link_to_Solution := new Solution'(s);

  begin
    if Is_Null(sols)
     then Construct(ls,sols);
     else temp := sols;
          while not Is_Null(temp) loop
            last := temp;
            temp := Tail_Of(temp);
          end loop;
          Construct(ls,tmp);
          Swap_Tail(last,tmp);
    end if;
  end Add;

  procedure Add ( sols : in out Solution_List; s : in Solution;
                  tol : in double_float; other : out natural ) is

    last,temp,tmp : Solution_List;
    ls : Link_to_Solution := new Solution'(s);
    s2 : Solution(s.n);
    count : natural := 1;

  begin
    other := 0;
    if Is_Null(sols)
     then Construct(ls,sols);
     else temp := sols;
          while not Is_Null(temp) loop
            s2 := Head_Of(temp).all;
            if Equal(s,s2,tol)
             then other := count;
                  Clear(ls);
                  return;
             else last := temp;
                  temp := Tail_Of(temp);
                  count := count + 1;
            end if;
          end loop;
          Construct(ls,tmp);
          Swap_Tail(last,tmp);
    end if;
  end Add;

  procedure Change ( sols : in out Solution_List; pos : in positive;
                     s : in Solution; tol : in double_float;
                     other : out natural ) is
  begin
    if pos <= Length_Of(sols)
     then declare
            temp : Solution_List := sols;
            ls : Link_to_Solution;
          begin
            other := 0;
            for i in 1..Length_Of(temp) loop
              ls := Head_Of(temp);
              if i = pos
               then ls.v := s.v;
                    ls.m := s.m;
                    ls.t := s.t;
                    Set_Head(temp,ls);
                    return;
               elsif Equal(s,ls.all,tol)
                   then other := i;
                        return;
              end if;
              temp := Tail_Of(temp);
            end loop;
          end;
    end if;
  end Change;

  procedure Set_Continuation_Parameter
               ( sols : in out Solution_List; t : in double_complex ) is

    tmp : Solution_List := sols;

  begin
    while not Is_Null(tmp) loop
      declare
        ls : Link_to_Solution := Head_Of(tmp);
      begin
        ls.t := t;
        Set_Head(tmp,ls);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Set_Continuation_Parameter;

  procedure Change_Multiplicity
                ( sols : in out Solution_List; pos : in positive;
                  m : in natural ) is
  begin
    if pos <= Length_Of(sols)
     then declare
            temp : Solution_List := sols;
            ls : Link_to_Solution;
          begin
            for i in 1..(pos-1) loop
              temp := Tail_Of(temp);
            end loop;
            ls := Head_Of(temp);
            ls.m := m;
            Set_Head(temp,ls);
          end;
    end if;
  end Change_Multiplicity;

  procedure Remove ( sols : in out Solution_List; pos : in positive ) is

    first,second,temp : Solution_List;
    ls : Link_to_Solution;

  begin
    if pos <= Length_Of(sols)
     then if pos = 1
           then if Is_Null(Tail_Of(sols))
                 then Clear(sols);
                 else ls := Head_Of(sols);
                      Clear(ls);
                      sols := Tail_Of(sols);
                end if;
           else second := sols;
                for i in 1..(pos-1) loop
                  first := second;
                  second := Tail_Of(first);
                end loop;
                ls := Head_Of(second);
                Clear(ls);
                temp := Tail_Of(second);
                Swap_Tail(first,temp);
          end if;
    end if;
  end Remove;

  procedure Delete ( sols : in out Solution_List ) is

    continue : boolean;

  begin
    continue := true;
    -- looking for the first element in sols that can stay :
    while not Is_Null(sols) and continue loop
      declare
        ls : Link_to_Solution := Head_Of(sols);
      begin
        if To_Be_Removed(ls.m)
         then Clear(ls);
              sols := Tail_Of(sols);
	 else continue := false;
        end if;
      end;
    end loop;
    if not Is_Null(sols)
     then -- first element of sols can stay in the list
	  declare
	    first,second : Solution_List;
          begin
	    first := sols;
	    second := Tail_Of(first);
	    while not Is_Null(second) loop
	      declare
		ls : Link_to_Solution := Head_Of(second);
		temp : Solution_List;
              begin
		if To_Be_Removed(ls.m)
		 then Clear(ls);
		      temp := Tail_Of(second);
		      Swap_Tail(first,temp);
                end if;
	      end;
	      first := second;
	      second := Tail_Of(first);
            end loop;
          end;
    end if;
  end Delete;
 
  procedure Remove_All ( sols : in out Solution_List; flag : in natural ) is

    continue : boolean;

  begin
    continue := true;
    -- looking for the first element in sols that can stay :
    while not Is_Null(sols) and continue loop
      declare
        ls : Link_to_Solution := Head_Of(sols);
      begin
        if ls.m = flag
         then Clear(ls);
              sols := Tail_Of(sols);
	 else continue := false;
        end if;
      end;
    end loop;
    if not Is_Null(sols)
     then -- first element of s can stay in the list
	  declare
	    first,second : Solution_List;
          begin
	    first := sols;
	    second := Tail_Of(first);
	    while not Is_Null(second) loop
	      declare
		ls : Link_to_Solution := Head_Of(second);
		temp : Solution_List;
              begin
		if ls.m = flag
		 then Clear(ls);
		      temp := Tail_Of(second);
		      Swap_Tail(first,temp);
                end if;
	      end;
	      first := second;
	      second := Tail_Of(first);
            end loop;
          end;
    end if;
  end Remove_All;
    
-- DESTRUCTORS :


  procedure Clear ( ls : in out Link_to_Solution ) is

    procedure free is new unchecked_deallocation(Solution,Link_to_Solution);

  begin
    free(ls);
  end Clear;

  procedure Shallow_Clear ( sl : in out Solution_List ) is
  begin
    List_of_Solutions.Clear(List_of_Solutions.List(sl));
  end Shallow_Clear;

  procedure Deep_Clear ( sl : in out Solution_List ) is

    temp : Solution_List := sl;
    ls : Link_to_Solution;

  begin
    while not Is_Null(temp) loop
      ls := Head_Of(temp);
      Clear(ls);
      temp := Tail_Of(temp);
    end loop;
    Shallow_Clear(sl);
  end Deep_Clear;

  procedure Clear ( sa : in out Solution_Array ) is
  begin
    for i in sa'range loop
      Clear(sa(i));
    end loop;
  end Clear;

end Solutions;
