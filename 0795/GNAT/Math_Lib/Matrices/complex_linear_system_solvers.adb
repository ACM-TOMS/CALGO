package body Complex_Linear_System_Solvers is

  use Complex_Vectors;

-- AUXLILIARIES :

  function cabs ( c : double_complex ) return double_float is
  begin
    return (ABS(REAL_PART(c)) + ABS(IMAG_PART(c)));
  end cabs;

  function dconjg ( x : double_complex ) return double_complex is
  begin
    return CMPLX(REAL_PART(x),-IMAG_PART(x));
  end dconjg;

  function csign ( x,y : double_complex ) return double_complex is
  begin
    return (CMPLX(cabs(x)) * y / CMPLX(cabs(y)));
  end csign;

-- TARGET ROUTINES :

  procedure Scale ( a : in out Matrix; b : in out Complex_Vectors.Vector ) is

    fac : double_complex;

    function Maximum ( a : in Matrix; i : in integer ) return double_complex is

      res : integer := a'first(2);
      max : double_float := cabs(a(i,res));
      tmp : double_float;

    begin
      for j in a'first(2)+1..a'last(2) loop
        tmp := cabs(a(i,j));
        if tmp > max
         then max := tmp; res := j;
        end if;
      end loop;
      return a(i,res);
    end Maximum;

    procedure Divide ( a : in out Matrix; b : in out Vector;
                       i : in integer; fac : in double_complex ) is
    begin
      for j in a'range(2) loop
        a(i,j) := a(i,j)/fac;
      end loop;
      b(i) := b(i)/fac;
    end Divide;
  
  begin
    for i in a'range(1) loop
      fac := Maximum(a,i);
      Divide(a,b,i,fac);
    end loop;
  end Scale;

-- TARGET ROUTINES :

  procedure lufac ( a : in out Matrix; n : in integer;
                    ipvt : out Integer_Vectors.Vector; info : out integer ) is

    kp1,l,nm1 : integer;
    smax : double_float;
    temp : double_complex;

  begin
    info := 0;
    nm1 := n - 1;
    if nm1 >= 1
     then for k in 1..nm1 loop
            kp1 := k + 1;

          -- find the pivot index l

            l := k; smax := cabs(a(k,k));  --modulus(a(k,k));
            for i in kp1..n loop
              if cabs(a(i,k)) > smax --modulus(a(i,k)) > smax
               then l := i;
                    smax := cabs(a(i,k)); --modulus(a(i,k));
              end if;
            end loop;
            ipvt(k) := l;

            if smax = 0.0
             then -- this column is already triangularized
                  info := k;
             else

                  -- interchange if necessary

                  if l /= k
                   then temp := a(l,k);
                        a(l,k) := a(k,k);
                        a(k,k) := temp;
                  end if;
 
                  -- compute multipliers

                  temp := -CMPLX(1.0)/a(k,k);
                  for i in kp1..n loop
                    a(i,k) := temp * a(i,k);
                  end loop;

                  -- row elimination with column indexing
  
                  for j in kp1..n loop
                    temp := a(l,j);
                    if l /= k
                     then a(l,j) := a(k,j);
                          a(k,j) := temp;
                    end if;
                    for i in kp1..n loop
                      a(i,j) := a(i,j) + temp * a(i,k);
                    end loop;
                  end loop;

            end if;
         end loop;
    end if;
    ipvt(n) := n;
    if modulus(a(n,n)) = 0.0
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

    z : Complex_Vectors.Vector(1..n);
    info,kb,kp1,l : integer;
    s,sm,sum,anorm,ynorm : double_float;
    ek,t,wk,wkm : double_complex;
    ipvtt : Integer_Vectors.Vector(1..n);

  begin
    anorm := 0.0;                                    -- compute 1-norm of a
    for j in 1..n loop
      sum := 0.0;
      for i in 1..n loop
        sum := sum + cabs(a(i,j));
      end loop;
      if sum > anorm
       then anorm := sum;
      end if; 
    end loop;
    lufac(a,n,ipvtt,info);                                        -- factor
    for i in 1..n loop
      ipvt(i) := ipvtt(i);
    end loop;
    ek := CMPLX(1.0);                              -- solve ctrans(u)*w = e
    for j in 1..n loop
      z(j) := CMPLX(0.0);
    end loop;
    for k in 1..n loop
      if cabs(z(k)) /= 0.0
       then ek := csign(ek,-z(k));
      end if;
      if cabs(ek-z(k)) > cabs(a(k,k)) 
       then s := cabs(a(k,k))/cabs(ek-z(k));
            z := CMPLX(s) * z;
            ek := CMPLX(s) * ek;
      end if;
      wk := ek - z(k);
      wkm := -ek - z(k);
      s := cabs(wk);
      sm := cabs(wkm);
      if cabs(a(k,k)) = 0.0 
       then wk := CMPLX(1.0);
            wkm := CMPLX(1.0);
       else wk := wk / dconjg(a(k,k));
            wkm := wkm / dconjg(a(k,k));
      end if;
      kp1 := k + 1;
       if kp1 <= n
        then for j in kp1..n loop
               sm := sm + cabs(z(j)+wkm*dconjg(a(k,j)));
               z(j) := z(j) + wk*dconjg(a(k,j));
               s := s + cabs(z(j));
             end loop;
             if s < sm
              then t := wkm - wk;
                   wk := wkm;
                   for j in kp1..n loop
                     z(j) := z(j) + t*dconjg(a(k,j));
                   end loop;
             end if;
       end if;
       z(k) := wk;
     end loop;
     sum := 0.0;
     for i in 1..n loop
       sum := sum + cabs(z(i));
     end loop;
     s := 1.0 / sum;
     z := CMPLX(s) * z;
     for k in 1..n loop                           -- solve ctrans(l)*y = w
       kb := n+1-k;
       if kb < n
        then t := CMPLX(0.0);
             for i in (kb+1)..n loop
               t := t + dconjg(a(i,kb))*z(i);
             end loop;
             z(kb) := z(kb) + t;
       end if;
       if cabs(z(kb)) > 1.0 
        then s := 1.0 / cabs(z(kb));
             z := CMPLX(s) * z;
       end if;
       l := ipvtt(kb);
       t := z(l);
       z(l) := z(kb);
       z(kb)  := t;
     end loop;
     sum := 0.0;
     for i in 1..n loop
       sum := sum + cabs(z(i));
     end loop;
     s := 1.0 / sum;
     z := CMPLX(s) * z;
     ynorm := 1.0;
     for k in 1..n loop                                    -- solve l*v = y
       l := ipvtt(k);
       t := z(l);
       z(l) := z(k);
       z(k) := t;
       if k < n
        then for i in (k+1)..n loop
               z(i) := z(i) + t * a(i,k);
             end loop;
       end if;
       if cabs(z(k)) > 1.0
        then s := 1.0 / cabs(z(k));
             z := CMPLX(s) * z;
             ynorm := s * ynorm;
       end if;
     end loop;
     sum := 0.0;
     for i in 1..n loop
       sum := sum + cabs(z(i));
     end loop;
     s := 1.0 / sum;
     z := CMPLX(s) * z;
     ynorm := s * ynorm;
     for k in 1..n loop                                    -- solve u*z = v
       kb := n+1-k;
       if cabs(z(kb)) > cabs(a(kb,kb))
        then s := cabs(a(kb,kb)) / cabs(z(kb));
             z := CMPLX(s) * z;
             ynorm := s * ynorm;
       end if;
       if cabs(a(kb,kb)) = 0.0
        then z(kb) := CMPLX(1.0);
        else z(kb) := z(kb) / a(kb,kb);
       end if;
       t := -z(kb);
       for i in 1..(kb-1) loop
         z(i) := z(i) + t * a(i,kb);
       end loop;
     end loop;
     sum := 0.0;                                       -- make znorm = 1.0
     for i in 1..n loop
       sum := sum + cabs(z(i));
     end loop;
     s := 1.0 / sum;
     z := CMPLX(s) * z;
     ynorm := s * ynorm;
     if anorm = 0.0
      then rcond := 0.0;
      else rcond := ynorm/anorm;
     end if;
  end lufco;

  procedure lusolve ( a : in Matrix; n : in integer;
                      ipvt : in Integer_Vectors.Vector;
                      b : in out Complex_Vectors.Vector ) is

    l,nm1,kb : integer;
    temp : double_complex;
 
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
    for k in 1..n loop                                     -- solve u*x = y
      kb := n+1-k;
      b(kb) := b(kb) / a(kb,kb);
      temp := -b(kb);
      for j in 1..(kb-1) loop
        b(j) := b(j) + temp * a(j,kb);
      end loop;
    end loop;
  end lusolve;

  procedure Triangulate ( a : in out Matrix; n,m : in integer ) is

    max,cbs : double_float;
    temp : double_complex;
    pivot,k,kcolumn : integer;
    tol : constant double_float := 10.0**(-10);

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := 0.0;                                             -- find pivot
      pivot := 0;
      for l in k..n loop
        cbs := cabs(a(l,kcolumn));
        if (cbs > tol) and then (cbs > max)
         then max := cbs;
              pivot := l;
        end if;
      end loop;
      if pivot = 0
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
            a(k,kcolumn) := CMPLX(1.0);
            for i in (k+1)..n loop
              for j in (kcolumn+1)..m loop
                a(i,j) := a(i,j) - a(i,kcolumn) * a(k,j);
              end loop;
              a(i,kcolumn) := CMPLX(0.0);
            end loop;
            k := k + 1;
            kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Triangulate;

  procedure Diagonalize ( a : in out Matrix; n,m : in integer ) is

    max : double_float;
    temp : double_complex;
    pivot,k,kcolumn : integer;

  begin
    k := 1;
    kcolumn := 1;
    while (k <= n) and (kcolumn <= m) loop
      max := 0.0;                                               -- find pivot
      for l in k..n loop
        if cabs(a(l,kcolumn)) > max
         then max := cabs(a(l,kcolumn));
              pivot := l;
        end if;
      end loop;
      if max = 0.0
       then kcolumn := kcolumn + 1;
       else if pivot /= k                        -- interchange if necessary
             then for i in 1..m loop
                    temp := a(pivot,i);
                    a(pivot,i) := a(k,i);
                    a(k,i) := temp;
                  end loop;
            end if;
            for j in (kcolumn+1)..m loop                    -- diagonalize a
              a(k,j) := a(k,j) / a(k,kcolumn);
            end loop;
            a(k,kcolumn) := CMPLX(1.0);
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
              a(j,kcolumn) := CMPLX(0.0);
            end loop;
            for j in (k+1)..n loop
              a(j,kcolumn) := CMPLX(0.0);
            end loop;
            k := k + 1;
            kcolumn := kcolumn + 1;
      end if;
    end loop;
  end Diagonalize;

end Complex_Linear_System_Solvers;
