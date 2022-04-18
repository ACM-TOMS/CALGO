with Greatest_Common_Divisor;  use Greatest_Common_Divisor;

package body Integer_Linear_System_Solvers is

-- SCALERS :

  function Divisors ( a : Matrix ) return Vector is

  -- DESCRIPTION :
  --   Returns a vector containing the gcd of the elements of each row.

    v : Vector(a'range(1));
   
  begin
    for i in a'range(1) loop
      v(i) := a(i,a'first(2));
      if v(i) /= 1
       then for j in (a'first(2)+1)..a'last(2) loop
	      v(i) := gcd(v(i),a(i,j));
	      exit when (v(i) = 1);
            end loop;
      end if;
    end loop;
    return v;
  end Divisors;

  function Scale ( a : Matrix ) return Matrix is

    v : Vector(a'range(1)) := Divisors(a);
    b : Matrix(a'range(1),a'range(2));

  begin
    for i in b'range(1) loop
      if (v(i) = 1) or (v(i) = 0)
       then for j in b'range(2) loop
	      b(i,j) := a(i,j);
            end loop;
       else for j in b'range(2) loop
	      b(i,j) := a(i,j) / v(i);
            end loop;
      end if;
    end loop;
    return b;
  end Scale;

  procedure Scale ( a : in out Matrix; v : out Vector ) is

    dv : Vector(a'range(1)) := Divisors(a);

  begin
    for i in a'range(1) loop
      if (dv(i) = 1) or (dv(i) = 0)
       then null;
       else for j in a'range(2) loop
	      a(i,j) := a(i,j) / dv(i);
            end loop;
      end if;
    end loop;
    v := dv;
  end Scale;

  procedure Scale ( a : in out Matrix ) is

    v : Vector(a'range(1)) := Divisors(a);

  begin
    for i in a'range(1) loop
      if (v(i) = 1) or (v(i) = 0)
       then null;
       else for j in a'range(2) loop
	      a(i,j) := a(i,j) / v(i);
            end loop;
      end if;
    end loop;
  end Scale;

  procedure Scale ( a : in out Matrix; row,col : in integer ) is

    g : integer := a(row,col);

  begin
    if g /= 1
     then for l in (col+1)..a'last(2) loop
            g := gcd(g,a(row,l));
            exit when g = 1;
          end loop;
    end if;
    if (g /= 0) and (g /= 1)
     then for l in row..a'last(2) loop
            a(row,l) := a(row,l)/g;
          end loop;
    end if;
  end Scale;

-- STATIC TRIANGULATORS :

  procedure Upper_Triangulate ( l : out Matrix; a : in out Matrix ) is

    row,pivot,temp,aa,bb,ka,lb,d,a_rowk,a_ik : integer;
    ll : Matrix(a'range(1),a'range(1));

  begin
    for i in ll'range(1) loop
      for j in ll'range(2) loop
        ll(i,j) := 0;
      end loop;
      ll(i,i) := 1;
    end loop;
    row := a'first(1);
    for j in a'first(2)..a'last(2) loop
      pivot := row-1;                                         -- find pivot
      for i in row..a'last(1) loop
	if a(i,j) /= 0
	 then pivot := i;
	      exit;
        end if;
      end loop;
      if pivot >= row 
       then if pivot /= row                      -- interchange if necessary
            then for k in a'range(2) loop
	           temp := a(row,k);
	           a(row,k) := a(pivot,k);
	           a(pivot,k) := temp;
	         end loop;
		 for k in ll'range(2) loop
		   temp := ll(row,k);
		   ll(row,k) := ll(pivot,k);
		   ll(pivot,k) := temp;
                 end loop;
           end if;
           for i in (row+1)..a'last(1) loop                  -- make zeroes
             if a(i,j) /= 0
              then aa := a(row,j); bb := a(i,j);
                   gcd(aa,bb,ka,lb,d);
                   aa := aa/d;     bb := bb/d;
                   if (aa = bb) and then ka = 0
                    then ka := lb; lb:= 0;
                   end if;
                   if (aa = -bb) and then ka = 0
                    then ka := -lb; lb := 0;
                   end if;
		   for k in j..a'last(2) loop
		     a_rowk := a(row,k); a_ik := a(i,k);
		     a(row,k) :=  ka*a_rowk + lb*a_ik;
		     a(i,k)   := -bb*a_rowk + aa*a_ik;
		   end loop;
		   for k in ll'range(2) loop
		     a_rowk := ll(row,k); a_ik := ll(i,k);
		     ll(row,k) :=  ka*a_rowk + lb*a_ik;
		     ll(i,k)   := -bb*a_rowk + aa*a_ik;
                   end loop;
             end if;
           end loop;
	   row := row + 1;
       end if;
       exit when row > a'last(1);
    end loop;
    l := ll;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( a : in out Matrix ) is

    row,pivot,temp,aa,bb,ka,lb,d,a_rowk,a_ik : integer;

  begin
    row := a'first(1);
    for j in a'first(2)..a'last(2) loop
      pivot := row-1;                                           -- find pivot
      for i in row..a'last(1) loop
	if a(i,j) /= 0
	 then pivot := i;
	      exit;
        end if;
      end loop;
      if pivot >= row 
       then if pivot /= row                       -- interchange if necessary 
            then for k in a'range(2) loop
	           temp := a(row,k);
	           a(row,k) := a(pivot,k);
	           a(pivot,k) := temp;
	         end loop;
           end if;
           for i in (row+1)..a'last(1) loop                    -- make zeroes
             if a(i,j) /= 0
              then aa := a(row,j); bb := a(i,j);
                   gcd(aa,bb,ka,lb,d);
                   aa := aa/d;     bb := bb/d;
                   if (aa = bb) and then ka = 0
                    then ka := lb; lb:= 0;
                   end if;
                   if (aa = -bb) and then ka = 0
                    then ka := -lb; lb := 0;
                   end if;
		   for k in j..a'last(2) loop
		     a_rowk := a(row,k);
		     a_ik   := a(i,k);
		     a(row,k) :=  ka*a_rowk + lb*a_ik;
		     a(i,k)   := -bb*a_rowk + aa*a_ik;
		   end loop;
             end if;
           end loop;
	   row := row + 1;
       end if;
       exit when row > a'last(1);
    end loop;
  end Upper_Triangulate;

  procedure Upper_Triangulate ( a : in out Matrix; ipvt : in out Vector ) is

    row,pivot,temp,aa,bb,ka,lb,d,a_rowk,a_ik : integer;

  begin
    row := a'first(1);
    for j in a'first(2)..a'last(2) loop
      pivot := row-1;                                           -- find pivot
      for i in row..a'last(1) loop
	if a(i,j) /= 0
	 then pivot := i;
	      exit;
        end if;
      end loop;
      if pivot >= row 
       then if pivot /= row                       -- interchange if necessary
            then for k in a'range(2) loop
	           temp := a(row,k);
	           a(row,k) := a(pivot,k);
	           a(pivot,k) := temp;
	         end loop;
                 temp := ipvt(row);
                 ipvt(row) := ipvt(pivot);
                 ipvt(pivot) := temp;
           end if;
           for i in (row+1)..a'last(1) loop                   -- make zeroes
             if a(i,j) /= 0
              then aa := a(row,j); bb := a(i,j);
                   gcd(aa,bb,ka,lb,d);
                   aa := aa/d;     bb := bb/d;
                   if (aa = bb) and then ka = 0
                    then ka := lb; lb:= 0;
                   end if;
                   if (aa = -bb) and then ka = 0
                    then ka := -lb; lb := 0;
                   end if;
		   for k in j..a'last(2) loop
		     a_rowk := a(row,k);
		     a_ik   := a(i,k);
		     a(row,k) :=  ka*a_rowk + lb*a_ik;
		     a(i,k)   := -bb*a_rowk + aa*a_ik;
		   end loop;
             end if;
           end loop;
	   row := row + 1;
       end if;
       exit when row > a'last(1);
    end loop;
  end Upper_Triangulate;

  procedure Lower_Triangulate ( a : in out Matrix; u : out Matrix ) is

    column,pivot,temp,aa,bb,ka,lb,d,a_kcolumn,a_kj : integer;
    uu : Matrix(a'range(2),a'range(2));

  begin

    for i in uu'range(1) loop
      for j in uu'range(2) loop
        uu(i,j) := 0;
      end loop;
      uu(i,i) := 1;
    end loop;
    column := a'first(2);
    for i in a'first(1)..a'last(1) loop
      pivot := column-1;                                      -- find pivot
      for j in column..a'last(2) loop
        if a(i,j) /= 0
         then pivot := j;
              exit;
        end if;
      end loop;
      if pivot >= column
       then if pivot /= column                  -- interchange if necessary
             then for k in a'range(1) loop
                    temp := a(k,column);
                    a(k,column) := a(k,pivot);
                    a(k,pivot) := temp;
                  end loop;
                  for k in uu'range(1) loop
                    temp := uu(k,column);
                    uu(k,column) := uu(k,pivot);
                    uu(k,pivot) := temp;
                  end loop;
            end if;
            for j in (column+1)..a'last(2) loop               -- make zeroes
              if a(i,j) /= 0
               then aa := a(i,column); bb := a(i,j);
                    gcd(aa,bb,ka,lb,d);
                    aa := aa/d;        bb := bb/d;
                    if (aa = bb) and then ka = 0
                     then ka := lb; lb:= 0;
                    end if;
                    if (aa = -bb) and then ka = 0
                     then ka := -lb; lb := 0;
                    end if;
                    for k in i..a'last(1) loop
                      a_kcolumn := a(k,column); a_kj := a(k,j);
                      a(k,column) := a_kcolumn*ka    + a_kj*lb;
                      a(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
                    end loop;
                    for k in uu'range(1) loop
                      a_kcolumn := uu(k,column); a_kj := uu(k,j);
                      uu(k,column) := a_kcolumn*ka    + a_kj*lb;
                      uu(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
                    end loop;
              end if;
            end loop;
            column := column + 1;
      end if;
      exit when column > a'last(2);
    end loop;
    u := uu;
  end Lower_Triangulate;

  procedure Lower_Triangulate ( a : in out Matrix ) is

    column,pivot,temp,aa,bb,ka,lb,d,a_kcolumn,a_kj : integer;

  begin
    column := a'first(2);
    for i in a'first(1)..a'last(1) loop
      pivot := column-1;                                      -- find pivot
      for j in column..a'last(2) loop
        if a(i,j) /= 0
         then pivot := j;
              exit;
        end if;
      end loop;
      if pivot >= column
       then if pivot /= column                  -- interchange if necessary
             then for k in a'range(1) loop
                    temp := a(k,column);
                    a(k,column) := a(k,pivot);
                    a(k,pivot) := temp;
                  end loop;
            end if;
            for j in (column+1)..a'last(2) loop               -- make zeroes
             if a(i,j) /= 0
              then aa := a(i,column); bb := a(i,j);
                   gcd(aa,bb,ka,lb,d);
                   aa := aa/d;        bb := bb/d;
                   if (aa = bb) and then ka = 0
                    then ka := lb; lb:= 0;
                   end if;
                   if (aa = -bb) and then ka = 0
                    then ka := -lb; lb := 0;
                   end if;
                   for k in i..a'last(1) loop
                     a_kcolumn := a(k,column);
                     a_kj := a(k,j);
                     a(k,column) := a_kcolumn*ka    + a_kj*lb;
                     a(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
                   end loop;
             end if;
           end loop;
           column := column + 1;
      end if;
      exit when column > a'last(2);
    end loop;
  end Lower_Triangulate;

  procedure Lower_Triangulate ( a : in out Matrix; ipvt : in out Vector ) is

    column,pivot,temp,aa,bb,ka,lb,d,a_kcolumn,a_kj : integer;

  begin
    column := a'first(2);
    for i in a'first(1)..a'last(1) loop
      pivot := column-1;                                       -- find pivot
      for j in column..a'last(2) loop
        if a(i,j) /= 0
         then pivot := j;
              exit;
        end if;
      end loop;
      if pivot >= column
       then if pivot /= column                   -- interchange if necessary
             then for k in a'range(1) loop
                    temp := a(k,column);
                    a(k,column) := a(k,pivot);
                    a(k,pivot) := temp;
                  end loop;
                  temp := ipvt(column);
                  ipvt(column) := ipvt(pivot);
                  ipvt(pivot) := temp;
            end if;
            for j in (column+1)..a'last(2) loop               -- make zeroes
              if a(i,j) /= 0
               then aa := a(i,column); bb := a(i,j);
                    gcd(aa,bb,ka,lb,d);
                    aa := aa/d;        bb := bb/d;
                    if (aa = bb) and then ka = 0
                     then ka := lb; lb:= 0;
                    end if;
                    if (aa = -bb) and then ka = 0
                     then ka := -lb; lb := 0;
                    end if;
                    for k in i..a'last(1) loop
                      a_kcolumn := a(k,column);
                      a_kj := a(k,j);
                      a(k,column) := a_kcolumn*ka    + a_kj*lb;
                      a(k,j)      := a_kcolumn*(-bb) + a_kj*aa;
                    end loop;
              end if;
            end loop;
            column := column + 1;
      end if;
      exit when column > a'last(2);
    end loop;
  end Lower_Triangulate;

-- SELECTORS :

  function Det ( a : Matrix ) return integer is

  -- NOTE :
  --   The triangulation is implemented independently to keep track
  --   of row interchanges.

    res : integer := 1;
    m : matrix(a'range(1),a'range(2));
    row,pivot,temp,aa,bb,ka,lb,d,m_rowk,m_ik : integer;

  begin
    m := a;   -- triangulate m
    row := m'first(1);
    for j in m'first(2)..m'last(2) loop
      pivot := row-1;                                           -- find pivot
      for i in row..m'last(1) loop
        if m(i,j) /= 0
         then pivot := i;
              exit;
        end if;
      end loop;
      if pivot >= row
       then if pivot /= row                       -- interchange if necessary
             then for k in m'range(2) loop
                    temp := m(row,k);
                    m(row,k) := m(pivot,k);
                    m(pivot,k) := temp;
                  end loop;
                  res := -res;
            end if;
            for i in (row+1)..m'last(1) loop                   -- make zeroes
              if m(i,j) /= 0
               then aa := m(row,j); bb := m(i,j);
                    gcd(aa,bb,ka,lb,d);
                    aa := aa/d;     bb := bb/d;
                    if (aa = bb) and then ka = 0
                     then ka := lb; lb:= 0;
                    end if;
                    if (aa = -bb) and then ka = 0
                     then ka := -lb; lb := 0;
                    end if;
                    for k in j..m'last(2) loop
                      m_rowk := m(row,k);
                      m_ik   := m(i,k);
                      m(row,k) :=  ka*m_rowk + lb*m_ik;
                      m(i,k)   := -bb*m_rowk + aa*m_ik;
                    end loop;
              end if;
            end loop;
            row := row + 1;
       end if;
       exit when row > m'last(1);
    end loop;
    for k in m'range(1) loop
      res := res*m(k,k);
    end loop;
    return res;
  end Det;

  function Per ( i,n : natural; a : matrix; kk : vector ) return natural is
  begin
    if i = n+1
     then return 1;
     else declare
            res : natural := 0;
            kkk : vector(kk'range) := kk;
          begin
            for j in kk'range loop
              if a(i,j) /= 0 and then kk(j) /= 0
               then kkk(j) := kkk(j) - 1;
                    res := res + a(i,j)*Per(i+1,n,a,kkk);
                    kkk(j) := kkk(j) + 1;
              end if;
            end loop;
            return res;
          end;
    end if;
  end Per;

  function Per ( i,n : natural; a : matrix; kk : vector; max : natural )
               return natural is
  begin
    if i = n+1
     then return 1;
     else declare
            res : natural := 0;
            kkk : vector(kk'range) := kk;
          begin
            for j in kk'range loop
              if a(i,j) /= 0 and then kk(j) /= 0
               then kkk(j) := kkk(j) - 1;
                    res := res + a(i,j)*Per(i+1,n,a,kkk,max);
                    kkk(j) := kkk(j) + 1;
              end if;
              exit when res >= max;
            end loop;
            return res;
          end;
    end if;
  end Per;

  function Per ( a : matrix; k : vector ) return natural is

  -- ALGORITHM :
  --   Row expansion without memory, as developed by C.W. Wampler,
  --   see `Bezout Number Calculations for Multi-Homogeneous Polynomial
  --        Systems', Appl. Math. Comput. 51:(2-3), 143-157, 1992.

  begin
    return Per(1,a'last(1),a,k);
  end Per;

  function Per ( a : matrix; k : vector; max : natural ) return natural is

  -- ALGORITHM :
  --   Row expansion without memory, as developed by C.W. Wampler,
  --   see `Bezout Number Calculations for Multi-Homogeneous Polynomial
  --        Systems', Appl. Math. Comput. 51:(2-3), 143-157, 1992.

  begin
    return Per(1,a'last(1),a,k,max);
  end Per;

  function Rank ( a : Matrix ) return natural is

    res : natural := 0; 
    m : Matrix(a'range(1),a'range(2));
    column : integer;

  begin
    m := a;
    Upper_Triangulate(m);
   -- compute the length of chain of nonzero elements in m :
   -- search first nonzero element in first row of m :
    column := m'first(2)-1;
    for k in m'range(2) loop
      if m(m'first(1),k) /= 0
       then column := k;
      end if;
      exit when (column = k);
    end loop;
    if column < m'first(2)
     then return 0;   -- all elements of m are zero
     else for k in m'range(1) loop
            exit when column > m'last(2);
            if m(k,column) /= 0
             then res := res + 1;
             else -- search for next nonzero element on row k :
                  for l in column+1..m'last(2) loop
                    if m(k,l) /= 0
                     then column := l;
                          res := res + 1;
                    end if;
                    exit when (column = l);
                  end loop;
            end if;
            column := column + 1;
          end loop;
    end if;
    return res;
  end Rank;

-- DYNAMIC TRIANGULATOR :

  procedure Triangulate ( l : in integer; m : in out matrix;
                          ipvt : in out vector; piv : out integer ) is

  -- DESCRIPTION :
  --   Updates lth row of m such that m remains upper triangular.

    pivot,tmp,a,b,lcmab,faca,facb : integer;
    index : integer;                -- first nonzero element in previous row
    tmpv : vector(m'range(2));

  begin
    Switch(ipvt,l,m);                      -- pivoting for previous unknowns
    index := 1;                         -- update : make l-1 zeroes in row l
    for k in 1..(l-1) loop
      if m(l,index) /= 0 and then m(k,index) /= 0    -- make m(l,index) zero
       then a := m(k,index); b := m(l,index);
            lcmab := lcm(a,b);
            if lcmab < 0 then lcmab := -lcmab; end if;
            facb := lcmab/b; faca := lcmab/a;
            if facb > 0
             then for i in index..m'last(2) loop
                    m(l,i) :=  facb*m(l,i) - faca*m(k,i);
                  end loop;
             else for i in index..m'last(2) loop
                    m(l,i) := -facb*m(l,i) + faca*m(k,i);
                  end loop;
            end if;
      end if;
      if m(k,index) /= 0
       then index := index + 1;
      end if;
    end loop;
    pivot := 0;                                              -- search pivot
    for k in l..m'last(2)-1 loop
      if m(l,k) /= 0
       then pivot := k;
      end if;
      exit when pivot /= 0;
    end loop;
    if pivot > l
     then for k in 1..l loop              -- interchange columns l and pivot
            tmp := m(k,l);
            m(k,l) := m(k,pivot);
            m(k,pivot) := tmp;
          end loop;
          tmp := ipvt(l);
          ipvt(l) := ipvt(pivot);
          ipvt(pivot) := tmp;
    end if;
    piv := pivot;
  end Triangulate;

  procedure Switch ( ipvt : in vector; index : in integer;
                     m : in out matrix ) is

    tmpv : Vector(m'range(2));

  begin
    for k in tmpv'range loop
      tmpv(k) := m(index,k);
    end loop;
    for k in tmpv'range loop
      m(index,k) := tmpv(ipvt(k));
    end loop;
  end Switch;

  procedure Switch ( ipvt : in vector; first,last : in integer;
                     m : in out matrix) is

    tmpv : vector(m'range(2));

  begin
    for index in first..last loop
      for k in tmpv'range loop
        tmpv(k) := m(index,k);
      end loop;
      for k in tmpv'range loop
        m(index,k) := tmpv(ipvt(k));
      end loop;
    end loop;
  end Switch;

  procedure Switch ( l,pivot,index : in integer; m : in out matrix ) is

    tmp : integer;

  begin
    if l /= pivot
     then tmp := m(index,l);
          m(index,l) := m(index,pivot);
          m(index,pivot) := tmp;
    end if;
  end Switch;

  procedure Switch ( l,pivot : in integer;
                     first,last : in integer; m : in out matrix ) is

    tmp : integer;

  begin
    if l /= pivot
     then for index in first..last loop
            tmp := m(index,l);
            m(index,l) := m(index,pivot);
            m(index,pivot) := tmp;
          end loop;
    end if;
  end Switch;

-- SOLVERS :

  function Check0 ( a : Matrix; x : Vector ) return boolean is
   
   -- DESCRIPTION :
   --   Returns true if x is a solution of the system a*x = 0.

    tmp : Vector(a'range(1));

  begin
    tmp := a*x;
    for i in tmp'range loop
      if tmp(i) /= 0
       then return false;
      end if;
    end loop;
    return true;
  end Check0;

  procedure Solve0 ( a : in Matrix; x : in out Vector ) is

  -- ALGORITHM :
  --   An intermediate, generating matrix tmp will be constructed,
  --   such that 
  --      1) the solution x to tmp*x = 0 is the same of a*x = 0;
  --      2) tmp(i,i) and tmp(i,ind) are the only nonzero entries.
  --   Before this construction, it will be checked whether there
  --   exists a zero column and the index ind must be determined.
  --   After the definition of tmp, the back substitution process
  --   yields a solution.

    piv,ind : integer;
    tmp : Matrix(a'first(1)..(a'last(1)+1),a'range(2));
    pivots,res : Vector(x'range);
    zero_column : boolean;

  begin
   -- initialization of tmp,ind and pivots :
    for i in tmp'range(1) loop
      for j in tmp'range(2) loop
	tmp(i,j) := 0;
      end loop;
    end loop;
    for i in pivots'range loop
      pivots(i) := i;
    end loop;
    ind := x'first(1)-1;
    for i in a'range(1) loop
      piv := pivots'first-1;
      for j in a'range(2) loop
	if a(i,j) /= 0
         then piv := pivots(j);
	      pivots(j) := pivots(i);
	      pivots(i) := piv;
	      exit;
        end if;
      end loop;
      zero_column := true;
      for j in a'first(1)..i loop
        tmp(j,i) := a(j,pivots(i));
	if zero_column and then tmp(j,i) /= 0
	 then zero_column := false;
        end if;
      end loop;
      if piv < pivots'first or else zero_column or else tmp(i,i) = 0
       then ind := i; exit;
      end if;
    end loop;
    if zero_column
     then x := (x'range => 0);
	  x(ind) := 1;
     elsif (ind < x'first(1)) and (a'last(1) >= a'last(2))
         then x := (x'range => 0);
         else
           if ind < x'first(1)
	    then ind := a'last(1)+1;
	         for j in tmp'range(2) loop
	           tmp(ind,j) := 0;
                 end loop;
	         zero_column := true;
	         for j in a'range(1) loop
	           tmp(j,ind) := a(j,pivots(ind));
	           if zero_column and then tmp(j,ind) /= 0
	            then zero_column := false;
                   end if;
                 end loop;
           end if;
           if zero_column
	    then x := (x'range => 0);
	         x(ind) := 1;
            else
             -- construct generating matrix :
              for i in reverse (tmp'first(2)+1)..(ind-1) loop  -- i = column
                for k in tmp'first(1)..(i-1) loop
                  if tmp(k,i) /= 0  -- make tmp(k,i) zero
                   then declare 
                          aa,bb,d : integer;
                        begin
                          aa := tmp(i,i); bb := tmp(k,i);
                          d := gcd(aa,bb);
                          aa := aa/d;     bb := bb/d;
                          for l in k..(i-1) loop
                            tmp(k,l) := aa*tmp(k,l);  -- tmp(i,l) = 0
                          end loop;
                          tmp(k,i) := 0; --aa*tmp(k,i) - bb*tmp(i,i);
                          tmp(k,ind) := aa*tmp(k,ind) - bb*tmp(i,ind);
                        end;
                        Scale(tmp,k,k);  -- to avoid numeric_error
                  end if;
                end loop; -- upper half of ith colum consists of zero entries
              end loop;
             -- generate x by back substitution :
              x(ind) := tmp(x'first,x'first);
              for i in (x'first+1)..(ind-1) loop
	        if tmp(i,i) /= 0
                 then x(ind) := lcm(tmp(i,i),x(ind));
                end if;
              end loop;
              for i in x'first..(ind-1) loop
	        if tmp(i,i) = 0
	         then x(i) := 0;
                 else x(i) := -(tmp(i,ind)*x(ind))/tmp(i,i);
                end if;
              end loop;
           end if;
    end if;
    res := (res'range => 0);                    -- take pivots into account
    for i in x'first..ind loop
      res(pivots(i)) := x(i);
    end loop;
    x := res;
  end Solve0;

  procedure Solve1 ( a : in Matrix; x : in out Vector;
	             b : in Vector; fail : out boolean ) is
  begin
    fail := false;
    for i in reverse x'range loop
      x(i) := b(i);
      for j in (i+1)..x'last loop
        x(i) := x(i) - a(i,j)*x(j);
      end loop;
      if x(i) /= 0 and then a(i,i) /= 0
       then if x(i) mod a(i,i) = 0
             then x(i) := x(i)/a(i,i);
             else fail := true; return;
            end if;
      end if;
    end loop;
  end Solve1;

  procedure Solve1 ( a : in Matrix; b : in out Vector;
		     fail : out boolean ) is
  begin
    fail := false;
    for i in reverse b'range loop
      for j in (i+1)..b'last loop
        b(i) := b(i) - a(i,j)*b(j);
      end loop;
      if b(i) /= 0 and then a(i,i) /= 0
       then if b(i) mod a(i,i) = 0
             then b(i) := b(i)/a(i,i);
             else fail := true; return;
            end if;
      end if;
    end loop;
  end Solve1;

  function Solve ( m : Matrix; ipvt : Vector ) return vector is

    x,res : vector(ipvt'range);
    a : matrix(m'first(1)..m'last(1)-1,m'range(2));
    ind : integer := a'first(1);       -- index for the current row number
    cnt0 : natural := 0;                 -- counts the number of zero rows

  begin
    for k in a'range(1) loop
      if m(k,k) /= 0                         -- otherwise : skip zero row !
       then for l in a'range(2) loop
              a(ind,l) := m(k,l);
            end loop;
            ind := ind + 1;
       else for l in a'range(2) loop
              a(a'last(1) - cnt0,l) := m(k,l);
            end loop;
            cnt0 := cnt0 + 1;
      end if;
    end loop;
    x := (x'range => 0);
    Solve0(a,x);
    for k in res'range loop
      res(ipvt(k)) := x(k);
    end loop;
    if res(res'last) < 0
     then return -res;
     else return res;
    end if;
  end Solve;

end Integer_Linear_System_Solvers;
