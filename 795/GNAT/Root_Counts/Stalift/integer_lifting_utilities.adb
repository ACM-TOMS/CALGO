with Integer_Vectors_Utilities;   use Integer_Vectors_Utilities;
with Power_Lists;                 use Power_Lists;

package body Integer_Lifting_Utilities is

  function Adaptive_Lifting ( l : Array_of_Lists ) return Vector is

    res : Vector(l'range);
    fac : constant natural := 3;     -- multiplication factor
    max : constant natural := 23;    -- upper bound for lifting

  begin
    for i in l'range loop
      res(i) := fac*Length_Of(l(i));
      if res(i) > max
       then res(i) := max;
      end if;
    end loop;
    return res;
  end Adaptive_Lifting;

  function Select_Subsystem ( p : Laur_Sys; mix : Vector; mic : Mixed_Cell )
                            return Laur_Sys is

    res : Laur_Sys(p'range);
    cnt : natural := 0;

  begin
    for k in mix'range loop
      for l in 1..mix(k) loop
        cnt := cnt + 1;
        res(cnt) := Select_Terms(p(cnt),mic.pts(k));
      end loop;
    end loop;
    return res;
  end Select_Subsystem;

  function Perform_Lifting ( n : natural; l : List; p : Poly ) return Poly is

    res : Poly := Null_Poly;
    tmp : List := l;

  begin
    while not Is_Null(tmp) loop
      declare
        d : Link_to_Vector := Head_Of(tmp);
        dr : Link_to_Vector := Reduce(d,n+1);
        t : Term;
      begin
        t.cf := Coeff(p,Degrees(dr));
        t.dg := Degrees(d);
        Plus_Term(res,t);
        Clear(dr);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Perform_Lifting;

  function Perform_Lifting
              ( n : natural; mix : Vector; lifted : Array_of_Lists;
                p : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(p'range);
    cnt : natural := 1;

  begin
    for k in mix'range loop
      for l in 1..mix(k) loop
        res(cnt) := Perform_Lifting(n,lifted(k),p(cnt));
        cnt := cnt+1;
      end loop;
    end loop;
    return res;
  end Perform_Lifting;

  function Copy_Lifting ( lifted : List; pt : Link_to_Vector )
                        return Link_to_Vector is

  -- DESCRIPTION :
  --   Searches the correspoinding point in the list lifted and returns
  --   the lifted point.  If the corresponding point has not been found,
  --   then the original point pt will be returned.

    tmp : List := lifted;
    lpt,res : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      lpt := Head_Of(tmp);
      if Equal(lpt(pt'range),pt.all)
       then res := new Integer_Vectors.Vector'(lpt.all);
            return res;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return pt;
  end Copy_Lifting;

  function Copy_Lifting ( lifted,pts : List ) return List is

  -- DESCRIPTION :
  --   Copies the lifting on the points lifted to the points in pts,
  --   i.e., each point in pts will get the same lifting as the corresponding
  --   lifted point in the list lifted.

    res : List;
    tmp : List := pts;

  begin
    while not Is_Null(tmp) loop
      Construct(Copy_Lifting(lifted,Head_Of(tmp)),res);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Copy_Lifting;

  procedure Search_Lifting ( l : in List; pt : in Vector;
                             found : out boolean; lif : out integer ) is

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
    res(res'last) := 1;
    res(res'last) := Conservative_Lifting(mixsub,k,res);
    return res;
  end Induced_Lifting;

  function Induced_Lifting
               ( n : natural; mix : Vector; points : Array_of_Lists;
                 mixsub : Mixed_Subdivision ) return Array_of_Lists is

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

  procedure Constant_Lifting
                ( l : in List; liftval : in natural;
                  lifted,lifted_last : in out List ) is

    tmp : List := l;
    pt : Link_to_Vector;
 
  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      declare
        lpt : Link_to_Vector := new Vector(pt'first..pt'last+1);
      begin
        lpt(pt'range) := pt.all;
        lpt(lpt'last) := liftval;
        Append(lifted,lifted_last,lpt);
      end;
      tmp := Tail_Of(tmp);
    end loop;
  end Constant_Lifting;

  procedure Constant_Lifting
               ( al : in Array_of_Lists; liftval : in natural;
                 lifted,lifted_last : in out Array_of_Lists ) is
  begin
    for i in al'range loop
      Constant_Lifting(al(i),liftval,lifted(i),lifted_last(i));
    end loop;
  end Constant_Lifting;

  function Conservative_Lifting
               ( mic : Mixed_Cell; k : natural; point : Vector )
               return integer is

    sp : integer := mic.nor*Head_Of(mic.pts(k));
    spp : integer := mic.nor.all*point;
    res : integer;

  begin
    if sp < spp
     then return point(point'last);
     else if mic.nor(mic.nor'last) = 0
           then res := point(point'last);
           else spp := spp - point(point'last)*mic.nor(mic.nor'last);
                res := (sp - spp)/mic.nor(mic.nor'last) + 1;
          end if;
          return res;
    end if;
  end Conservative_Lifting;

  function Conservative_Lifting ( mixsub : Mixed_Subdivision; k : natural;
                                  point : Vector ) return integer is

    tmp : Mixed_Subdivision := mixsub;
    pt : Vector(point'range) := point;
    res : integer;

  begin
    while not Is_Null(tmp) loop
      pt(pt'last) := Conservative_Lifting(Head_Of(tmp),k,pt);
      tmp := Tail_Of(tmp);
    end loop;
    res := pt(pt'last);
    Clear(pt);
    return res;
  end Conservative_Lifting;

  function Lower_Lifting ( mic : Mixed_Cell; k : natural; point : Vector )
                         return integer is
  begin
    if Is_In(mic.pts(k),point)
     then return 0;
     else declare
            pt : Vector(point'range) := point;
          begin
            pt(pt'last) := 0;
            return Conservative_Lifting(mic,k,pt);
          end;
    end if;
  end Lower_Lifting;

  function Lower_Lifting ( mixsub : Mixed_Subdivision; k : natural;
                           point : Vector ) return integer is

    lif : integer := point(point'last);
    tmp : Mixed_Subdivision := mixsub;
    max : integer := 0;

  begin
    while not Is_Null(tmp) loop
      lif := Lower_Lifting(Head_Of(tmp),k,point);
      if lif > max
       then max := lif;
      end if;
      exit when max = point(point'last);
      tmp := Tail_Of(tmp);
    end loop;
    return max;
  end Lower_Lifting;

end Integer_Lifting_Utilities;
