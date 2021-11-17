with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Float_Matrices;                     use Float_Matrices;
with Dictionaries;
with Linear_Programming;                 use Linear_Programming;

with Integer_Vectors_Utilities;          use Integer_Vectors_Utilities;
with Integer_Vectors_of_Vectors;
with Lists_of_Vectors_Utilities;         use Lists_of_Vectors_Utilities;
with Transformations;                    use Transformations;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Face_Enumerators;                   use Face_Enumerators;

--with text_io,integer_io;                 use text_io,integer_io;
--with Integer_Vectors_io;                 use Integer_Vectors_io;
--with Float_Vectors_io;                   use Float_Vectors_io;
--with Lists_of_Integer_Vectors_io;        use Lists_of_Integer_Vectors_io;

package body Vertices is

  function Is_In_Hull ( point : Vector; l : List ) return boolean is

    fpt : Float_Vectors.Vector(point'range);

  begin
    for i in point'range loop
      fpt(i) := double_float(point(i));
    end loop;
    return Is_In_Hull(fpt,l);
  end Is_In_Hull;

  function Is_In_Hull ( point : Float_Vectors.Vector; l : List )
                      return boolean is

  -- ALGORITHM:
  --   The following linear program will be solved:
  --
  --    min u_0 + u_1 + .. + u_n
  --
  --        l_1*p_1i + l_2*p_2i + .. + l_m*p_mi + u_i*q_i = q_i   i=1,2,..,n
  --
  --        l_1      + l_2      + .. + l_m      + u_0     = 1
  --
  --  to determine whether q belongs to the convex hull spanned by the
  --  vectors p_1,p_2,..,p_m.
  --  If all u_i are zero and all constraints are satisfied,
  --  then q belongs to the convex hull.

  -- CONSTANTS :

    n  : constant natural := point'last;           
    m  : constant natural := 2*n+2;               -- number of constraints
    nb : constant natural := Length_Of(l)+n+1;    -- number of unknowns
    eps : constant double_float := 10.0**(-10);

  -- VARIABLES :

    dic : matrix(0..m,0..nb);
    sol : Float_Vectors.vector(1..nb);
    inbas : Integer_Vectors.Vector(1..m) := (1..m => 0);
    outbas : Integer_Vectors.Vector(1..nb) := (1..nb => 0);
    nit : natural := 0;
    feasi : boolean;

    tmp : List;
    pt : Link_to_Vector;
    s : double_float;

  begin

   -- INITIALIZATION OF target :

    for i in 0..(nb-n-1) loop
      dic(0,i) :=  0.0;             -- sum of the lambda's
    end loop;
    for i in (nb-n)..nb loop
      dic(0,i) := -1.0;             -- sum of the mu's
    end loop;

   -- INITIALIZATION OF dic :

    for i in 0..(nb-n) loop
      dic(m-1,i) :=  1.0;           -- sum of the lambda's + mu_0
      dic(m,i)   := -1.0;
    end loop;
    for i in (nb-n+1)..nb loop
      dic(m-1,i) := 0.0;
      dic(m,i)   := 0.0;
    end loop;
    for i in 1..n loop
      for j in (nb-n)..nb loop
        if i /= j-nb+n
         then dic(i,j)   := 0.0;
              dic(i+n,j) := 0.0;
        end if;
      end loop;
    end loop;

    tmp := l;
    for j in 1..(nb-n-1) loop
      pt := Head_Of(tmp);
      for i in pt'range loop
        dic(i,j)   :=  double_float(pt(i));
        dic(i+n,j) := -double_float(pt(i));
      end loop;
      tmp := Tail_Of(tmp);
    end loop;

    for i in point'range loop
      dic(i,0)   :=  point(i);
      dic(i+n,0) := -point(i);
      dic(i,i+nb-n)   :=  point(i);
      dic(i+n,i+nb-n) := -point(i);
    end loop;

  -- SOLVE THE LINEAR PROGRAM :

    Dictionaries.Init_Basis(inbas,outbas);
    Dual_Simplex(dic,eps,inbas,outbas,nit,feasi);
    if not feasi
     then return false;
     else
       sol := Dictionaries.Primal_Solution(dic,inbas,outbas);

  -- CHECK THE SOLUTION :

       s := 0.0;
       for i in 1..(nb-n-1) loop
         s := s + sol(i);
       end loop;
       if abs(s - 1.0) > eps
        then return false;
       end if;
       s := 0.0;
       for i in (nb-n)..nb loop
         s := s + sol(i);
       end loop;
       if abs(s) > eps
        then return false;
       end if;

       return true;

    end if;

  end Is_In_Hull;

  function Vertex_Points ( l : List ) return List is

    result,result_last : List;
    tmp,rest,origrest,rest_last : List;
    pt : Link_to_Vector;

  begin
    if Is_Null(l) or else Is_Null(Tail_Of(l))
     then return l;
     else Copy(Tail_Of(l),rest);
          origrest := rest;
          rest_last := rest;
          while not Is_Null(Tail_Of(rest_last)) loop
            rest_last := Tail_Of(rest_last);
          end loop;
          tmp := l;
          for i in 1..Length_Of(l) loop
            pt := Head_Of(tmp);
            if not Is_In_Hull(pt.all,rest)
             then Append(result,result_last,pt.all);
                  Append(rest,rest_last,pt.all);
            end if;
            rest := Tail_Of(rest);
            tmp := Tail_Of(tmp);
          end loop;
          Clear(origrest);
          return result;
    end if;
  end Vertex_Points;

  function Vertex_Points1 ( l : List ) return List is

    len : constant natural := Length_Of(l);
    pts : Integer_Vectors_of_Vectors.Vector(1..len) := Shallow_Create(l);
    result,result_last : List;

    procedure Collect_Vertex ( i : in integer; cont : out boolean ) is
    begin
      Append(result,result_last,pts(i).all);
      cont := true;
    end Collect_Vertex;
    procedure Enum_Vertices is new Enumerate_Vertices(Collect_Vertex);

  begin
    Enum_Vertices(pts);
    return result;
  end Vertex_Points1;

  procedure Add ( pt : Link_to_Vector; l : in out List ) is

   -- DESCRIPTION :
   --   Constructs the point to the list, but without sharing.

    newpt : Link_to_Vector := new Vector'(pt.all);

  begin
    Construct(newpt,l);
  end Add;

  function Extremal_Points ( l : List; v : Link_to_Vector ) return List is

    min,max,sp : integer;
    tmp,res : List;
    minpt,maxpt,pt : Link_to_Vector;

  begin
    if Length_Of(l) <= 1
     then Copy(l,res);
     else pt := Head_Of(l);
          min := pt*v; max := min;
          minpt := pt; maxpt := pt;
          tmp := Tail_Of(l);
          while not Is_Null(tmp) loop
            pt := Head_Of(tmp);
            sp := pt*v;
            if sp > max
             then max := sp; maxpt := pt;
             elsif sp < min
                 then min := sp; minpt := pt;
            end if;
            tmp := Tail_Of(tmp);
          end loop;
          Add(minpt,res);
          if min /= max
           then Add(maxpt,res);
          end if;
    end if;
    return res;
  end Extremal_Points;

  function Extremal_Points ( k,n : natural; exl,l : List ) return List is

    res,tmp,nres : List;
    iv : Integer_Vectors_of_Vectors.Vector(1..k);
    v : Link_to_Vector := new Vector(1..n);
    sp : integer;
    pt : Link_to_Vector;
    done : boolean;

  begin
    tmp := exl;
    for j in iv'range loop
      iv(j) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    v := Compute_Normal(iv);
    sp := Head_Of(exl)*v;
    nres := Extremal_Points(l,v);
    tmp := nres;
    res := exl;
    done := false;
    while not done and not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt*v /= sp
       then Add(pt,res);
            done := true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    Clear(v);
    return res;
  end Extremal_Points;

  function Max_Extremal_Points ( k,n : natural; exl,l : List ) return List is

    res,tmp,nres : List;
    iv : Integer_Vectors_of_Vectors.Vector(1..k);
    v : Link_to_Vector := new Vector(1..n);
    sp : integer;
    pt : Link_to_Vector;
    done : boolean;

  begin
    tmp := exl;
    for j in iv'range loop
      iv(j) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    v := Compute_Normal(iv);
    sp := Head_Of(exl)*v;
    nres := Extremal_Points(l,v);
    --put("The computed normal : "); put(v); new_line;
    --put("The inner product : "); put(sp,1); new_line;
    --put_line("The list of extremal points :"); put(nres);
    tmp := nres;
    Copy(exl,res);
    done := false;
    while not done and not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt*v /= sp
       then Add(pt,res);
            done := true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    Clear(nres);
    if not done
     then -- for all points x in l : <x,v> = sp
          if n >= 2
           then declare
                  i : integer := Pivot(v);
                  t,invt : Transfo;
                  texl,tl,tres : List;
                begin
                  if i <= v'last
                   then t := Build_Transfo(v,i);
                        invt := Invert(t);
                        texl := Transform_and_Reduce(t,i,exl);
                        tl := Transform_and_Reduce(t,i,l);
                        tres := Max_Extremal_Points(k,n-1,texl,tl);
                      --put("The normal : "); put(v); new_line;
                      --put_line("The list of points : "); put(l);
                      --put_line("The list of extremal points :"); put(exl);
                      --put_line("The transformed list of points : "); put(tl);
                      --put_line("The transformed list of extremal points : "); 
                      --put(texl);
                      --put_line("The computed transformed extremal points : ");
                      --put(tres);
                      --put("k : "); put(k,1); new_line;
                        Clear(t); Clear(texl); Clear(tl);
                        if Length_Of(tres) = k+1
                         then res := Insert_and_Transform(tres,i,sp,invt);
                         else Copy(exl,res);
                        end if;
                      --put_line("The computed extremal points : "); put(res);
                        Clear(tres); --responsible for segmentation violation...
                        Clear(invt);
                   else Copy(exl,res);
                  end if;
                end;
           else Copy(exl,res);
          end if;
    end if;
    Clear(v);
    --put_line("The new list of extremal points : "); put(res);
    return res;
  end Max_Extremal_Points;

  function Extremal_Points ( n : natural; l : List ) return List is

    res : List;
    nor : Link_to_Vector := new Vector'(1..n => 1);
    k : natural;
  
  begin
    res := Extremal_Points(l,nor);
    Clear(nor);
    if Length_Of(res) < 2
     then return res;
     else k := 2;
          while k < n+1 loop
            res := Extremal_Points(k,n,res,l);
            exit when Length_Of(res) = k;
            k := k+1;
          end loop;
          return res;
    end if;
  end Extremal_Points;

  function Two_Extremal_Points ( n : natural; l : List ) return List is

   -- DESCRIPTION :
   --   Computes two extremal points of the list l.

    res,tmp : List;
    pt,minpt,maxpt : Link_to_Vector;
    min,max : integer;

  begin
    if Length_Of(l) <= 2
     then Copy(l,res);
     else for i in 1..n loop
            pt := Head_Of(l);
            max := pt(i); min := max;
            maxpt := pt;  minpt := pt;
            tmp := Tail_Of(l);
            while not Is_Null(tmp) loop
              pt := Head_Of(tmp);
              if pt(i) < min
               then min := pt(i); minpt := pt;
               elsif pt(i) > max
                   then max := pt(i); maxpt := pt;
              end if;
              tmp := Tail_Of(tmp);
            end loop;
            if Is_Null(res)
             then Add(minpt,res);
                  if min /= max
                   then Add(maxpt,res);
                  end if;
             else pt := Head_Of(res);
                  if pt.all /= minpt.all
                   then Add(minpt,res);
                   elsif pt.all /= maxpt.all
                       then Add(maxpt,res);
                  end if; 
            end if;
            exit when (Length_Of(res) = 2);
          end loop;
    end if;
    return res;
  end Two_Extremal_Points;

  function Max_Extremal_Points ( n : natural; l : List ) return List is

    res,wl,wres,newres : List;
    k : natural;
    nullvec,shiftvec : Link_to_Vector;
   -- shifted : boolean;

  begin
    if Length_Of(l) <= 2
     then Copy(l,res);
     else res := Two_Extremal_Points(n,l);
          k := Length_Of(res);
          if k = 2
           then --nullvec := new Vector'(1..n => 0);
                --if Is_In(nullvec,res)
                -- then Copy(l,wl);
                --      Copy(res,wres);
                --      shifted := false;
                -- else shiftvec := Head_Of(res);
                --      wl := Shift(shiftvec,l);
                --      wres := Shift(shiftvec,res);
                --      shifted := true;
                --end if;
                --Clear(nullvec);
                Copy(res,wres);
                while k < n+1 loop
                  newres := Max_Extremal_Points(k,n,wres,l);
                  Copy(newres,wres);
                  exit when Length_Of(wres) = k;
                  k := k+1;
                end loop;
                res := wres;
               -- Clear(res);
               -- if shifted
               --  then Min_Vector(shiftvec);
               --       res := Shift(shiftvec,wres);
               --       Clear(wres);
               --  else Copy(wres,res);
               -- end if;
               -- Clear(wl);
          end if;
    end if;
    return res;
  end Max_Extremal_Points;

end Vertices;
