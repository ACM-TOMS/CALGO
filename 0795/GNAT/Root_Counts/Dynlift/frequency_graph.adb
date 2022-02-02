package body Frequency_Graph is

-- CREATORS :

  function Occurrences ( i : natural; l : List ) return natural is

    res : natural := 0;
    tmp : List := l;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt(i) /= 0
       then res := res + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Occurrences;

  function Graph ( n : natural; supports : Array_of_Lists ) return matrix is

    res : matrix(1..n,supports'range);

  begin
    for i in 1..n loop
      for j in supports'range loop
        res(i,j) := Occurrences(i,supports(j));
      end loop;
    end loop;
    return res;
  end Graph;

-- MODIFIER :

  procedure Ignore ( m : in out matrix; point : in vector ) is
  begin
    for i in point'range loop
      if point(i) /= 0
       then for j in m'range(2) loop
              m(i,j) := 1;
            end loop;
      end if;
    end loop;
  end Ignore;

-- SELECTORS :

  function Occurrence ( i : natural; m : matrix ) return natural is

    res : natural := 0;

  begin
    for j in m'range(2) loop
      if m(i,j) /= 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Occurrence;

  function Occurrence ( i : natural; m : matrix; col : natural; perm : vector )
                      return natural is

    res : natural := 0;

  begin
    for j in col+1..m'last(2) loop
      if m(i,perm(j)) /= 0
       then res := res + 1;
      end if;
    end loop;
    return res;
  end Occurrence;

  function Occurrence ( v : vector; m : matrix ) return natural is

    min : natural := 1000000;
    occ : natural;

  begin
    for i in v'range loop
      if v(i) /= 0
       then occ := Occurrence(i,m);
            if occ < min
             then min := occ;
            end if;
      end if;
    end loop;
    return min;
  end Occurrence;

  function Occurrence ( v : vector; m : matrix; col : natural; perm : vector )
                      return natural is

    min : natural := 1000000;
    occ : natural;

  begin
    for i in v'range loop
      if v(i) /= 0
       then occ := Occurrence(i,m,col,perm);
            if occ < min
             then min := occ;
            end if;
      end if;
    end loop;
    return min;
  end Occurrence;

  function Lowest_Occurrence
               ( vec : Integer_Vectors_of_Vectors.Vector;
                 start : natural; m : matrix ) return natural is

    res : natural := start;
    min : natural := Occurrence(vec(start).all,m);
    occ : natural;

  begin
    for i in start+1..vec'last loop
      occ := Occurrence(vec(i).all,m);
      if occ < min
       then min := occ; res := i;
      end if;
    end loop;
    return res;
  end Lowest_Occurrence;

  function Lowest_Occurrence
               ( vec : Integer_Vectors_of_Vectors.Vector;
                 start : natural; m : matrix;
                 col : natural; perm : vector ) return natural is

    res : natural := start;
    min : natural := Occurrence(vec(start).all,m,col,perm);
    occ : natural;

  begin
    for i in start+1..vec'last loop
      occ := Occurrence(vec(i).all,m,col,perm);
      if occ < min
       then min := occ; res := i;
      end if;
    end loop;
    return res;
  end Lowest_Occurrence;

-- CONSTRUCTORS :

  function Sort ( l : List; m : matrix ) return List is

    res : List;
    vec : Integer_Vectors_of_Vectors.Vector(1..Length_Of(l))
           := Deep_Create(l);
    tmp : Link_to_Vector;
    low : natural;

  begin
    if Length_Of(l) <= 1
     then Copy(l,res);
     else for i in vec'first..vec'last-1 loop
            low := Lowest_Occurrence(vec,i,m);
            if low /= i
             then tmp := vec(i);
                  vec(i) := vec(low);
                  vec(low) := tmp;
            end if;
          end loop;
          res := Deep_Create(vec);
          Integer_Vectors_of_Vectors.Clear(vec);
    end if;
    return res;
  end Sort;

  procedure Sort ( l : in out List; m : in matrix ) is

    res : List := Sort(l,m);

  begin
    Copy(res,l); Deep_Clear(res);
  end Sort;

  function Sort ( l : List; m : matrix; col : natural; perm : vector )
                return List is

    res : List;
    vec : Integer_Vectors_of_Vectors.Vector(1..Length_Of(l))
           := Deep_Create(l);
    tmp : Link_to_Vector;
    low : natural;

  begin
    if Length_Of(l) <= 1
     then Copy(l,res);
     else for i in vec'first..vec'last-1 loop
            low := Lowest_Occurrence(vec,i,m,col,perm);
            if low /= i
             then tmp := vec(i);
                  vec(i) := vec(low);
                  vec(low) := tmp;
            end if;
          end loop;
          res := Deep_Create(vec);
          Integer_Vectors_of_Vectors.Clear(vec);
    end if;
    return res;
  end Sort;

  procedure Sort ( l : in out List; m : in matrix;
                   col : in natural; perm : in vector ) is

    res : List := Sort(l,m,col,perm);

  begin
    Copy(res,l); Deep_Clear(res);
  end Sort;

end Frequency_Graph;
