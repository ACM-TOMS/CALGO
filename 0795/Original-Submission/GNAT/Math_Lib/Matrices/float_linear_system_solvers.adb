package body Float_Linear_System_Solvers is

  use Float_Vectors;

-- AUXILIARY :

  function csign (x,y : double_float) return double_float is
  begin
    return abs(x) * y / abs(y);
  end csign;

-- STATIC TRIANGULATORS :

  procedure lufac ( a : in out matrix; n : in integer;
                    ipvt : out Integer_Vectors.Vector; info : out integer ) is

    kp1,l,nm1 : integer;
    smax : double_float;
    temp : double_float;

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1
     then for k in 1..nm1 loop
            kp1 := k + 1;
            l := k; smax := abs(a(k,k));            -- find the pivot index l
            for i in kp1..n loop
              if abs(a(i,k)) > smax 
               then l := i;
                    smax := abs(a(i,k)); 
              end if;
            end loop;
            ipvt(k) := l;
            if smax = 0.0
             then info := k;           -- this column is already triangulated
             else if l /= k                       -- interchange if necessary
                   then temp := a(l,k);
                        a(l,k) := a(k,k);
                        a(k,k) := temp;
                  end if;
                  temp := -1.0/a(k,k);                 -- compute multipliers
                  for i in kp1..n loop
                    a(i,k) := temp*a(i,k);
                  end loop;
                  for j in kp1..n loop                     -- row elimination
                    temp := a(l,j);
                    if l /= k
                     then a(l,j) := a(k,j);
                          a(k,j) := temp;
                    end if;
                    for i in kp1..n loop
                      a(i,j) := a(i,j) + temp*a(i,k);
                    end loop;
                  end loop;
            end if;
          end loop;
    end if;
    ipvt(n) := n;
    if a(n,n) = 0.0
     then info := n;
    end if;
  end lufac;

  procedure lufco ( a : in out Matrix; n : in integer;
                    ipvt : out Integer_Vectors.Vector;
                    rcond : out double_float ) is

  -- NOTE :
  --   rcond = 1/(norm(a)*(estimate of norm(inverse(a))))
  --   estimate = norm(z)/norm(y) where a*z = y and ctrans(a)*y = e.
  --   ctrans(a) is the conjugate transpose of a.
  --   The components of e are chosen to cause maximum local
  --   growth in teh elements of w where ctrans(u)*w = e.
  --   The vectors are frequently rescaled to avoid overflow.

    z : Float_Vectors.Vector(1..n);
    info,kb,kp1,l : integer;
    s,sm,sum,anorm,ynorm : double_float;
    ek,t,wk,wkm : double_float;
    ipvtt : Integer_Vectors.Vector(1..n);

  begin
    anorm := 0.0;                                     -- compute 1-norm of a
    for j in 1..n loop
      sum := 0.0;
      for i in 1..n loop
        sum := sum + abs(a(i,j));
      end loop;
      if sum > anorm
       then anorm := sum;
      end if; 
    end loop;
    lufac(a,n,ipvtt,info);                                         -- factor
    for i in 1..n loop
      ipvt(i) := ipvtt(i);
    end loop;
    ek := 1.0;                                     -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := 0.0;
    end loop;
    for k in 1..n loop
      if abs(z(k)) /= 0.0
       then ek := csign(ek,-z(k));
      end if;
      if abs(ek-z(k)) > abs(a(k,k)) 
       then s := abs(a(k,k))/abs(ek-z(k));
            z := s*z;
            ek := s * ek;
      end if;
      wk := ek - z(k);
      wkm := -ek - z(k);
      s := abs(wk);
      sm := abs(wkm);
      if abs(a(k,k)) = 0.0 
       then wk := 1.0;
            wkm := 1.0;
       else wk := wk / a(k,k);
            wkm := wkm / a(k,k);
      end if;
      kp1 := k + 1;
      if kp1 <= n
       then for j in kp1..n loop
              sm := sm + abs(z(j)+wkm*a(k,j));
              z(j) := z(j) + wk*a(k,j);
              s := s + abs(z(j));
            end loop;
            if s < sm
             then t := wkm - wk;
                  wk := wkm;
                  for j in kp1..n loop
                    z(j) := z(j) + t*a(k,j);
                  end loop;
            end if;
      end if;
      z(k) := wk;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + abs(z(i));
    end loop;
    s := 1.0 / sum;
    z := s*z;
    for k in 1..n loop                              -- solve ctrans(l)*y = w
      kb := n+1-k;
      if kb < n
       then t := 0.0;
            for i in (kb+1)..n loop
              t := t + a(i,kb)*z(i);
            end loop;
            z(kb) := z(kb) + t;
      end if;
      if abs(z(kb)) > 1.0 
       then s := 1.0 / abs(z(kb));
            z := s*z;
      end if;
      l := ipvtt(kb);
      t := z(l);
      z(l) := z(kb);
      z(kb)  := t;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + abs(z(i));
    end loop;
    s := 1.0 / sum;
    z := s*z;
    ynorm := 1.0;
    for k in 1..n loop                                      -- solve l*v = y
      l := ipvtt(k);
      t := z(l);
      z(l) := z(k);
      z(k) := t;
      if k < n
       then for i in (k+1)..n loop
              z(i) := z(i) + t * a(i,k);
            end loop;
      end if;
      if  abs(z(k)) > 1.0
       then s := 1.0 / abs(z(k));
            z := s*z;
            ynorm := s * ynorm;
      end if;
    end loop;
    sum := 0.0;
    for i in 1..n loop
      sum := sum + abs(z(i));
    end loop;
    s := 1.0 / sum;
    z := s*z;
    ynorm := s * ynorm;
    for k in 1..n loop                                      -- solve u*z = v
      kb := n+1-k;
      if abs(z(kb)) > abs(a(kb,kb))
       then s := abs(a(kb,kb)) / abs(z(kb));
            z := s*z;
            ynorm := s * ynorm;
      end if;
      if abs(a(kb,kb)) = 0.0
       then z(kb) := 1.0;
       else z(kb) := z(kb) / a(kb,kb);
      end if;
      t := -z(kb);
      for i in 1..(kb-1) loop
        z(i) := z(i) + t * a(i,kb);
      end loop;
    end loop;
    sum := 0.0;                                         -- make znorm = 1.0
    for i in 1..n loop
      sum := sum + abs(z(i));
    end loop;
    s := 1.0 / sum;
    z := s*z;
    ynorm := s * ynorm;
    if anorm = 0.0
     then rcond := 0.0;
     else rcond := ynorm/anorm;
    end if;
  end lufco;

  procedure lusolve ( a : in matrix; n : in integer;
                      ipvt : in Integer_Vectors.Vector;
                      b : in out Float_Vectors.vector ) is

    l,nm1,kb : integer;
    temp : double_float;
 
  begin
    nm1 := n-1;
    if nm1 >= 1                                             -- solve l*y = b
     then for k in 1..nm1 loop
            l := ipvt(k);
            temp := b(l);
            if l /= k 
             then b(l) := b(k);
                  b(k) := temp;
            end if;
            for i in (k+1)..n loop
              b(i) := b(i) + temp * a(i,k);
            end loop;
          end loop;
    end if;
    for k in 1..n loop                                      -- solve u*x = y
      kb := n+1-k;
      b(kb) := b(kb) / a(kb,kb);
      temp := -b(kb);
      for j in 1..(kb-1) loop
        b(j) := b(j) + temp * a(j,kb);
      end loop;
    end loop;
  end lusolve;

  procedure Triangulate ( a : in out Matrix; n,m : in integer ) is

    max : double_float;
    temp : double_float;
    pivot,k,kcolumn : integer;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := 0.0;                                              -- find pivot
      for l in k..n loop
        if abs(a(l,kcolumn)) > max
         then max := abs(a(l,kcolumn));
              pivot := l;
        end if;
      end loop;
      if max = 0.0
       then kcolumn := kcolumn + 1;
       else if pivot /= k                       -- interchange if necessary
             then for i in 1..m loop
                    temp := a(pivot,i);
                    a(pivot,i) := a(k,i);
                    a(k,i) := temp;
                  end loop;
            end if;
            for j in (kcolumn+1)..m loop                   -- triangulate a
              a(k,j) := a(k,j) / a(k,kcolumn);
            end loop;
            a(k,kcolumn) := 1.0;
            for i in (k+1)..n loop
              for j in (kcolumn+1)..m loop
                a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
              end loop;
            end loop;
            for j in (k+1)..n loop
              a(j,kcolumn) := 0.0;
            end loop;
            k := k + 1;
            kcolumn := kcolumn + 1;
       end if;
    end loop;
  end Triangulate;

  procedure Diagonalize ( a : in out Matrix; n,m : in integer ) is

    max,temp : double_float;
    pivot,k,kcolumn : integer;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := 0.0;                                              -- find pivot
      for l in k..n loop
        if abs(a(l,kcolumn)) > max
         then max := abs(a(l,kcolumn));
              pivot := l;
        end if;
      end loop;
      if max = 0.0
       then kcolumn := kcolumn + 1;
       else if pivot /= k                       -- interchange if necessary
             then for i in 1..m loop
                    temp := a(pivot,i);
                    a(pivot,i) := a(k,i);
                    a(k,i) := temp;
                  end loop;
            end if;
            for j in (kcolumn+1)..m loop                   -- diagonalize a
              a(k,j) := a(k,j) / a(k,kcolumn);
            end loop;
            a(k,kcolumn) := 1.0;
            for i in 1..(k-1) loop
              for j in (kcolumn+1)..m loop
                a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
              end loop;
            end loop;
            for i in (k+1)..n loop
              for j in (kcolumn+1)..m loop
                a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
              end loop;
            end loop;
            for j in 1..(k-1) loop
              a(j,kcolumn) := 0.0;
            end loop;
            for j in (k+1)..n loop
              a(j,kcolumn) := 0.0;
            end loop;
            k := k + 1;
            kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Diagonalize;

-- DYNAMIC TRIANGULATORS :

  procedure Upper_Triangulate
               ( row : in natural; mat : in out Matrix; tol : in double_float;
                 ipvt : in out Integer_Vectors.Vector; pivot : out integer ) is

   factor,tmp,max : double_float;
   piv,tpi : integer := 0;

  begin
    for j in mat'first(1)..row-1 loop
      if abs(mat(row,j)) > tol                 -- make mat(row,j) zero
       then factor := mat(row,j)/mat(j,j);
            for k in j..mat'last(2) loop
              mat(row,k) := mat(row,k) - factor*mat(j,k);
            end loop;
      end if;
    end loop;
    for j in row..ipvt'last loop             -- search pivot
      tmp := abs(mat(row,j));
      if tmp > tol
       then if piv = 0
             then max := tmp; piv := j;
             elsif tmp > max
                 then max := tmp; piv := j;
            end if;
      end if;
    end loop;
    pivot := piv;
    if piv /= 0                                -- zero row
     then if piv /= row                        -- interchange columns
           then for k in mat'first(1)..row loop
                  tmp := mat(k,row); mat(k,row) := mat(k,piv);
                  mat(k,piv) := tmp;
                end loop;
                tpi := ipvt(row); ipvt(row) := ipvt(piv); ipvt(piv) := tpi;
          end if;
    end if;
  end Upper_Triangulate;

  procedure Upper_Triangulate
               ( roweli : in natural; elim : in Matrix; tol : in double_float;
                 rowmat : in natural; mat : in out Matrix ) is

    factor : double_float;

  begin
    if abs(mat(rowmat,roweli)) > tol
     then factor := mat(rowmat,roweli)/elim(roweli,roweli);
          for i in roweli..mat'last(2) loop
            mat(rowmat,i) := mat(rowmat,i) - factor*elim(roweli,i);
          end loop;
    end if;
  end Upper_Triangulate;

  procedure Upper_Triangulate
               ( roweli : in natural; elim : in Matrix; tol : in double_float;
                 firstrow,lastrow : in natural; mat : in out Matrix ) is
  begin
    for i in firstrow..lastrow loop
      Upper_Triangulate(roweli,elim,tol,i,mat);
    end loop;
  end Upper_Triangulate;

  procedure Switch ( ipvt : in Integer_Vectors.Vector; row : in integer;
                     mat : in out Matrix ) is

    tmp : Float_Vectors.Vector(mat'range(2));

  begin
    for k in tmp'range loop
      tmp(k) := mat(row,k);
    end loop;
    for k in ipvt'range loop
      mat(row,k) := tmp(ipvt(k));
    end loop;
    for k in ipvt'last+1..mat'last(2) loop
      mat(row,k) := tmp(k);
    end loop;
  end Switch;

  procedure Switch ( k,pivot,first,last : in integer; mat : in out Matrix ) is

    tmp : double_float;

  begin
    if k /= pivot
     then for i in first..last loop
            tmp := mat(i,k);
            mat(i,k) := mat(i,pivot);
            mat(i,pivot) := tmp;
          end loop;
    end if;
  end Switch;

  function Solve ( mat : Matrix; tol : double_float;
                   ipvt : Integer_Vectors.Vector )
                 return Float_Vectors.Vector is

    res,x : Float_Vectors.Vector(mat'range(2)) := (mat'range(2) => 0.0);
    index : integer;

  begin
    for i in mat'range(1) loop
      index := i;
      exit when i > mat'last(2);
      exit when abs(mat(i,i)) < tol;
    end loop;
    if (abs(mat(index,index)) > tol) and then (index < mat'last(2))
     then index := index + 1;
    end if;
    x(index) := 1.0;
    for i in reverse mat'first(1)..(index-1) loop
      x(i) := -mat(i,index);
      for j in i+1..index-1 loop
        x(i) := x(i) - mat(i,j)*x(j);
      end loop;
      x(i) := x(i)/mat(i,i);
    end loop;
    for k in ipvt'range loop
      res(ipvt(k)) := x(k);
    end loop;
    for k in ipvt'last+1..res'last loop
      res(k) := x(k);
    end loop;
    return res;
  end Solve;

end Float_Linear_System_Solvers;
