package body Degrees_in_Sets_of_Unknowns is

  function Degree ( t : Term; s : Set ) return integer is

    sum : integer := 0;

  begin
    if Extent_Of(s) > 0
     then for i in t.dg'range loop
            if Is_In(s,i)
             then sum := sum + t.dg(i);
            end if;
          end loop;
    end if;
    return sum;
  end Degree;

  function Degree ( p : Poly; s : Set ) return integer is

    res : integer := -1;

    procedure Degree_Term ( t : in Term; continue : out boolean ) is
      sum : integer := Degree(t,s);
    begin
      if sum > res
       then res := sum;
      end if;
      continue := true;
    end Degree_Term;
    procedure Degree_Terms is new Visiting_Iterator(Degree_Term);

  begin
    Degree_Terms(p);
    return res;
  end Degree;

  function Degree_Table ( p : Poly_Sys; z : Partition ) return matrix is

    res : matrix(p'range,z'range);

  begin
    for i in p'range loop
      for j in z'range loop
        res(i,j) := Degree(p(i),z(j));
      end loop;
    end loop;
    return res;
  end Degree_Table;

end Degrees_in_Sets_of_Unknowns;
