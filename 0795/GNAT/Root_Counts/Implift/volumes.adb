--with text_io,Integer_Vectors_io;       use text_io,Integer_Vectors_io;

with Integer_Vectors;                    use Integer_Vectors;
with Integer_Vectors_Utilities;          use Integer_Vectors_Utilities;
with Integer_Vectors_of_Vectors;
with Integer_Matrices;                   use Integer_Matrices;
with Integer_Linear_System_Solvers;      use Integer_Linear_System_Solvers;
with Transformations;                    use Transformations;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Integer_Support_Functions;          use Integer_Support_Functions;

with Lists_of_Vectors_Utilities;         use Lists_of_Vectors_Utilities;
with Arrays_of_Lists_Utilities;          use Arrays_of_Lists_Utilities;

with Face_Enumerators;
with Face_Enumerator_of_Sum;

package body Volumes is

-- INVARIANT CONDITION :
--   In order to get p(v) > 0, the zero vector must be in the first list,
--   so shifting is necessary.
--   The shift vector must be equal to the maximal element in the list,
--   w.r.t. the graded lexicographic ordening.
--   In this way, the shift vector is unique, which allows to `mirror'
--   the same operations for the mixed homotopy continuation.

-- AUXILIARIES :

  function Create_Facet 
             ( n : natural; facet : Vector;
               pts : Integer_Vectors_of_Vectors.Vector )
             return Integer_Vectors_of_Vectors.Vector is

    res : Integer_Vectors_of_Vectors.Vector(1..n);
    cnt : natural := 0;

  begin
    for i in facet'range loop
      cnt := cnt+1;
      res(cnt) := pts(facet(i));
    end loop;
    return res;
  end Create_Facet;

  function Determinant ( vv : Integer_Vectors_of_Vectors.Vector )
                       return integer is

    a : matrix(vv'range,vv'range);

  begin
    for k in a'range(1) loop
      for l in a'range(2) loop
        a(k,l) := vv(k)(l);
      end loop;
    end loop;
   -- Upper_Triangulate(a);
    return Det(a);
  end Determinant;

-- VOLUME COMPUTATIONS :

  function Vol ( n : natural; l : List ) return natural is

  -- DESCRIPTION :
  --   Computes the volume of the simplex spanned by the list
  --   and the origin.

    res : integer;
    vv : Integer_Vectors_of_Vectors.Vector(1..n);
    tmp : List;

  begin
    tmp := l;
    for i in vv'range loop
      vv(i) := Head_Of(tmp);
      tmp := Tail_Of(tmp);
    end loop;
    res := Determinant(vv);
    if res < 0
     then return -res;
     else return res;
    end if;
  end Vol;

  function Volume ( n,i : natural; l : List; v : Vector; pv : integer )
                  return natural is

  -- DESCRIPTION :
  --   The points in l all belong to the same hyperplane:
  --    < x , v > - pv = 0;
  --   this procedure computes the volume of the polytope generated by
  --   the points in l and the origin.

    ll : natural := Length_Of(l);

  begin
    if ll < n
     then return 0;
     elsif ll = n
	 then return Vol(n,l);
         else declare
	        t : Transfo; 
 	        tl,rl : List;
		vl : integer;
              begin
	        if pv > 0
 	         then t := Build_Transfo(v,i);
                      tl := t*l;
	              rl := Reduce(tl,i);
                      Clear(t); Deep_Clear(tl);
                     -- Remove_Duplicates(rl);
                      if Length_Of(rl) >= n-1
                       then vl := pv*Volume(n-1,rl);
                       else vl := 0;
                      end if;
		      Deep_Clear(rl);
		      return vl;
                 else return 0;
                end if;
              end;
    end if;
  end Volume;

  function Sum ( n,m : natural; l : List ) return natural is

  -- DESCRIPTION :
  --   Computes the volume of the polytope generated by the points in l;
  --   where n > 1 and n < m = Length_Of(l).

    res : natural;
    fix : Link_to_Vector;
    pts : Integer_Vectors_of_Vectors.Vector(1..m-1);
    nulvec : Vector(1..n) := (1..n => 0);
    vl,wl,vl_last : List;
    sh : boolean;

    procedure Compute_Facet ( facet : in vector; cont : out boolean ) is

      vv : constant Integer_Vectors_of_Vectors.Vector
              := Create_Facet(n,facet,pts);
      pl : List;
      v : Link_to_Vector;
      i,pv,sup : integer;

    begin
      v := Compute_Normal(vv);
      i := Pivot(v);
      if i <= v'last
       then pv := vv(vv'first)*v;
	    if pv < 0
	     then for j in v'range loop
		    v(j) := -v(j);
                  end loop;
		  pv := -pv;
            end if;
	    if (pv > 0) and then not Is_In(vl,v)
             then sup := Maximal_Support(wl,v.all);
		  if sup = pv
		   then Append(vl,vl_last,v);
			pl := Face(wl,v.all,pv);
			res := res + Volume(n,i,pl,v.all,pv);
			Deep_Clear(pl);
                  end if;
            end if;
      end if;
      cont := true;
    end Compute_Facet;

    procedure Compute_Facets is
      new Face_Enumerators.Enumerate_Faces(Compute_Facet);
   
  begin
    res := 0;
    if not Is_In(l,nulvec)
     then fix := Graded_Max(l);
          wl := Shift(l,fix); sh := true;
     else wl := l;            sh := false;
    end if;
    Move_to_Front(wl,nulvec);
    pts := Shallow_Create(Tail_Of(wl));
    Compute_Facets(n-1,pts);
    Deep_Clear(vl);
    if sh
     then Deep_Clear(wl); Clear(fix);
    end if;
    return res;
  end Sum;

  function Volume ( n : natural; l : List ) return natural is
    m : natural;
  begin
    if n = 1
     then declare
            min,max : integer := 0;
	    d : Link_to_Vector := Head_Of(l);
          begin
	    Min_Max(l,d'first,min,max);
	    return (max - min);
          end;
     else m := Length_Of(l);
	  if m <= n
	   then return 0;
	   else return Sum(n,m,l);
          end if;
    end if;
  end Volume;

  function Volume ( n,i : natural; l : List; v : Vector;
                    pv : integer; tv : Tree_of_Vectors ) return natural is

  -- DESCRIPTION :
  --   The points in l all belong to the same hyperplane:
  --    < x , v > - pv = 0;
  --   this procedure computes the volume of the polytope generated by
  --   the points in l and the origin.

    ll : natural := Length_Of(l);

  begin
    if ll < n
     then return 0;
     elsif ll = n
	 then return Vol(n,l);
         else declare
	        t : Transfo; 
 	        tl,rl : List;
		vl : integer;
              begin
	        if pv > 0
 	         then t := Build_Transfo(v,i);
                      tl := t*l;
	              rl := Reduce(tl,i);
                      Clear(t); Deep_Clear(tl);
                     -- Remove_Duplicates(rl);
                      if Length_Of(rl) >= n-1
                       then vl := pv*Volume(n-1,rl,tv);
                       else vl := 0;
                      end if;
		      Deep_Clear(rl);
		      return vl;
                 else return 0;
                end if;
              end;
    end if;
  end Volume;

  function Sum ( n,m : natural; l : List; tv : Tree_of_Vectors )
	       return natural is

  -- DESCRIPTION :
  --   Computes the volume of the polytope generated by the points in l;
  --   where n > 1 and n < m = Length_Of(l).
  --   The tree of degrees tv is not empty.

    res : natural;
    fix : Link_to_Vector;
    nulvec : Vector(1..n) := (1..n => 0);
    wl : List;
    tmp : Tree_of_Vectors;
    sh : boolean;

  begin
    res := 0;
    if not Is_In(l,nulvec)
     then fix := Graded_Max(l);
          wl := Shift(l,fix); sh := true;
     else wl := l;            sh := false;
    end if;
    Move_to_Front(wl,nulvec);
    tmp := tv;
    while not Is_Null(tmp) loop
      declare
	nd : node := Head_Of(tmp);
	v : Link_to_Vector := nd.d;
	pv,i : integer;
	pl : List;
      begin
	i := Pivot(v);
	pv := Maximal_Support(wl,v.all);
	pl := Face(wl,v.all,pv);
	if (nd.ltv = null) or else Is_Null(nd.ltv.all)
	 then res := res + Volume(n,i,pl,v.all,pv);
	 else res := res + Volume(n,i,pl,v.all,pv,nd.ltv.all);
        end if;
	Deep_Clear(pl);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    if sh
     then Deep_Clear(wl); Clear(fix);
    end if;
    return res;
  end Sum;

  function Volume ( n : natural; l : List; tv : Tree_of_Vectors ) 
                  return natural is
    m : natural;
  begin
    if n = 1
     then declare
            min,max : integer := 0;
	    d : Link_to_Vector := Head_Of(l);
          begin
	    Min_Max(l,d'first,min,max);
	    return (max - min);
          end;
     else m := Length_Of(l);
	  if m <= n
	   then return 0;
	   elsif not Is_Null(tv)
	       then return Sum(n,m,l,tv);
	       else return Sum(n,m,l);
          end if;
    end if;
  end Volume;

  procedure Volume ( n,i : in natural; l : in List; v : in Vector;
                     pv : in integer; tv : in out Tree_of_Vectors;
                     vlm : out natural ) is

  -- DESCRIPTION :
  --   The points in l all belong to the same hyperplane:
  --    < x , v > - pv = 0;
  --   this procedure computes the volume of the polytope generated by
  --   the points in l and the origin.

    ll : natural := Length_Of(l);
    vl : natural;

  begin
    if ll < n
     then vlm := 0;
     elsif ll = n
	 then vl := Vol(n,l);
              if vl > 0
	       then declare
                      nd : node;
                    begin
		      nd.d := new Integer_Vectors.Vector'(v);
		      Construct(nd,tv);
                    end;
              end if;
	      vlm := vl;
         else declare
	        t : Transfo; 
 	        tl,rl : List;
              begin
	        if pv > 0
 	         then t := Build_Transfo(v,i);
                      tl := t*l;
                      rl := Reduce(tl,i);
                      Clear(t); Deep_Clear(tl);
                     -- Remove_Duplicates(rl);
                      if Length_Of(rl) >= n-1
                       then declare
                              tmp : Tree_of_Vectors;
                            begin
                              Volume(n-1,rl,tmp,vl);
                              if vl = 0
                               then Clear(tmp);
                               else
                                 declare
                                   nd : node;
                                 begin
                                   nd.d := new Integer_Vectors.Vector'(v);
                                   if not Is_Null(tmp)
                                    then nd.ltv := new Tree_of_Vectors'(tmp);
                                   end if;
                                   Construct(nd,tv);
                                 end;
                              end if;
                            end;
                            vlm := pv*vl;
                       else vlm := 0;
                      end if;
                      Deep_Clear(rl);
                 else vlm := 0;
                end if;
              end;
    end if;
  end Volume;

  procedure Sum ( n,m : in natural; l : in List;
		  tv : in out Tree_of_Vectors; vlm : out natural ) is

  -- DESCRIPTION :
  --   Computes the volume of the polytope generated by the points in l;
  --   where n > 1 and n < m = Length_Of(l).

    res : natural;
    pts : Integer_Vectors_of_Vectors.Vector(1..m-1);
    fix : Link_to_Vector;
    nulvec : Vector(1..n) := (1..n => 0);
    wl : List;
    sh : boolean;

    procedure Compute_Facet ( facet : in vector; cont : out boolean ) is

      vv : constant Integer_Vectors_of_Vectors.Vector
                              := Create_Facet(n,facet,pts);
      pl : List;
      v : Link_to_Vector;
      i,pv,sup,pvol : integer;

    begin
      v := Compute_Normal(vv);
      i := Pivot(v);
      if i <= v'last
       then pv := vv(vv'first)*v;
	    if pv < 0
	     then for j in v'range loop
		    v(j) := -v(j);
                  end loop;
		  pv := -pv;
            end if;
	    if (pv > 0) and then not Is_In(tv,v)
             then sup := Maximal_Support(wl,v.all);
		  if sup = pv
		   then pl := Face(wl,v.all,pv);
			Volume(n,i,pl,v.all,pv,tv,pvol);
			res := res + pvol;
			Deep_Clear(pl);
                  end if;
            end if;
      end if;
      cont := true;
    end Compute_Facet;

    procedure Compute_Facets is
      new Face_Enumerators.Enumerate_Faces(Compute_Facet);
   
  begin
    res := 0;
    if not Is_In(l,nulvec)
     then fix := Graded_Max(l);
          wl := Shift(l,fix); sh := true;
     else wl := l;            sh := false;
    end if;
    Move_to_Front(wl,nulvec);
    pts := Shallow_Create(Tail_Of(wl));
    Compute_Facets(n-1,pts);
    if sh
     then Deep_Clear(wl); Clear(fix);
    end if;
    vlm := res;
  end Sum;

  procedure Volume ( n : in natural; l : in List; 
		     tv : in out Tree_of_Vectors; vlm : out natural ) is
    m : natural;
  begin
    if n = 1
     then declare
            min,max : integer := 0;
	    d : Link_to_Vector := Head_Of(l);
          begin
	    Min_Max(l,d'first,min,max);
	    vlm := max - min;
          end;
     else m := Length_Of(l);
	  if m <= n
	   then vlm := 0;
	   else Sum(n,m,l,tv,vlm);
          end if;
    end if;
  end Volume;

-- MIXED VOLUME COMPUTATIONS WITHOUT TREE OF USEFUL DIRECTIONS :

  function Two_Terms_Mixed_Vol ( n : natural; al : Array_of_Lists )
			       return natural is

  -- DESCRIPTION :
  --   returns the mixed volume of the polytopes generated by the
  --   points in al, where the first polytope is generated by only
  --   two points.

    first,second : Link_to_Vector;
    l : List renames al(al'first);
    res : natural;

  begin
    first := Head_Of(l);
    second := Head_Of(Tail_Of(l));
    declare
      d : Vector(first'range);
      piv : integer := 0;
    begin
      for i in d'range loop
	d(i) := first(i) - second(i);
	if (piv = 0) and then (d(i) /= 0)
	 then piv := i;
        end if;
      end loop;
      if piv = 0
       then return 0;
       else if d(piv) < 0
	     then Min_Vector(d);
            end if;
	    declare
	      t : Transfo := Rotate(d,piv);
	      tr_al : Array_of_Lists(al'first..(al'last-1));
              degen : boolean := false;
            begin
              Apply(t,d);
	      for i in tr_al'range loop
		tr_al(i) := Transform_and_Reduce(t,piv,al(i+1));
                Remove_Duplicates(tr_al(i));
                degen := (Length_Of(tr_al(i)) <= 1);
                exit when degen;
              end loop;
              Clear(t);
              if not degen
	       then res := d(piv)*Mixed_Volume(n-1,tr_al);
               else res := 0;
              end if;
              Deep_Clear(tr_al);
            end;
      end if;
    end;
    return res;
  end Two_Terms_Mixed_Vol;

  function Facet_Normal
               ( n : natural; facet,pts : Integer_Vectors_of_Vectors.Vector )
               return Vector is

    res,pt1,pt2 : Vector(1..n);
    im : matrix(1..n,1..n);
    cnt : natural := 0;

  begin
    for i in facet'range loop
      if facet(i)'length > 1
       then pt1 := pts(facet(i)(facet(i)'first)).all;
            for j in facet(i)'first+1..facet(i)'last loop
              pt2 := pts(facet(i)(j)).all;
              cnt := cnt + 1;
              for k in pt1'range loop
                im(cnt,k) := pt2(k) - pt1(k);
              end loop;
            end loop;
      end if;
    end loop;
    for j in 1..n loop
      im(n,j) := 0;
    end loop;
    Upper_Triangulate(im);
    Scale(im);
    res := (1..n => 0);
    Solve0(im,res);
    Normalize(res);
   -- put("The facet normal : "); put(res); new_line;
    return res;
  end Facet_Normal;
     
  function Minkowski_Sum ( n : natural; al : Array_of_Lists )
                         return natural is

  -- DESCRIPTION :
  --   Computes the mixed volume of the polytope generated
  --   by the points in al, where n > 1.

    res,m,ptslen : natural;
    vl,vl_last,al1 : List;
    typ,ind : Vector(1..n);
    perm,mix : Link_to_Vector;
    wal : Array_of_Lists(al'range) := Interchange2(al);

    procedure Update ( v : in Vector; i : in integer;
		       added : in out boolean ) is

    -- DESCRIPTION :
    --   This procedure computes the support of the first list
    --   in the direction v; if this support is not zero,
    --   the projection will be computed.
    --   The parameter added becomes true if v has been added to vl.

      pal : Array_of_Lists(al'first..(al'last-1));
      pv : integer;
      degen : boolean;

    begin
      if not Is_In(vl,v)
       then pv := Maximal_Support(al1,v);
            if pv > 0
             then Projection(wal,v,i,pal,degen);
	          if not degen
	           then declare
                          mv : integer := Mixed_Volume(n-1,pal);
                        begin
                          if mv > 0
   	                   then res := res + pv*mv;
                                Append(vl,vl_last,v);
                                added := true;
                          end if;
                        end;
                        Deep_Clear(pal);
                  end if;
            end if;
      end if;
    end Update;

    procedure Enumerate_Facets ( lpts : in Array_of_Lists; len : in natural ) is

      pts : Integer_Vectors_of_Vectors.Vector(1..len);
      cnt : integer;

      procedure Compute_Facet
                   ( facet : in Integer_Vectors_of_Vectors.Vector;
                     cont : out boolean ) is

        v,w : Vector(1..n);
        i : integer;
        added : boolean;

      begin
        v := Facet_Normal(n,facet,pts);
        i := Pivot(v);
        if i <= v'last
         then added := false;
              Update(v,i,added);
              if not added
	       then Min_Vector(v); w := v;
	       else w := -v; added := false;
              end if;
	      Update(w,i,added);
        end if;
        cont := true;
      end Compute_Facet;

      procedure Compute_Facets is 
        new Face_Enumerator_of_Sum.Enumerate_Faces_of_Sum(Compute_Facet);
       -- new Face_Enumerators.Enumerate_Faces_of_Sum(Compute_Facet);

    begin
      pts(ind(1)..ind(2)-1) := Shallow_Create(Tail_Of(al1));
      cnt := lpts'first + mix(mix'first);
      for i in mix'first+1..mix'last loop
        if i < mix'last
         then pts(ind(i)..ind(i+1)-1) := Shallow_Create(lpts(cnt));
              cnt := cnt + mix(i);
         else pts(ind(i)..len) := Shallow_Create(lpts(cnt));
        end if;
      end loop;
      Compute_Facets(ind(mix'range),typ(mix'range),n-1,pts);
    end Enumerate_Facets;
   
  begin
    m := Length_Of(wal(wal'first));
    if m = 2
     then return Two_Terms_Mixed_Vol(n,wal);
     elsif m > 2
         then
           Mixture(al,perm,mix);
          -- put("Type of mixture : "); put(mix); new_line;
          -- put(" with permutation : "); put(perm); new_line;
           wal := Permute(perm.all,al);
           declare
	     shiftvec : Link_to_Vector;
             nulvec : Vector(1..n) := (1..n => 0);
	     shifted : boolean;
             cnt : integer;
	   begin
	    -- SHIFT OF THE FIRST LIST ( then all pv >= 0) :
             if not Is_In(wal(wal'first),nulvec)
              then shiftvec := Graded_Max(wal(wal'first));
	           al1 := Shift(wal(wal'first),shiftvec); shifted := true;
              else al1 := wal(wal'first);                 shifted := false;
             end if;
            -- ENUMERATE FACES OF SUM OF THE RIGHT TYPE :
             Move_to_Front(al1,nulvec);
             wal(wal'first) := al1;
	     res := 0;
             typ(1) := mix(mix'first)-1;
             typ(2) := mix(mix'first+1); ind(1) := 1;
             ind(2) := ind(1) + Length_Of(al1) - 1;  -- skip null vector
             cnt := wal'first + mix(mix'first);
	     for i in mix'first+2..mix'last loop
               typ(i) := mix(i);
               ind(i) := ind(i-1) + Length_Of(wal(cnt));
               cnt := cnt + mix(i-1);
             end loop;
             ptslen := ind(mix'last) + Length_Of(wal(cnt)) - 1;
             Enumerate_Facets(wal,ptslen);
            -- CLEANING UP :
             Deep_Clear(vl); Clear(perm); Clear(mix);
             if shifted 
	      then Clear(shiftvec);
                   Deep_Clear(al1);
             end if;
             return res;
           end;
      else -- m < 2
           return 0;
    end if;
  end Minkowski_Sum;

  function Mixed_Volume ( n : natural; al : Array_of_Lists )
			return natural is
  begin
    if (n = 0) or Is_Null(al(al'first))
     then return 0;
     elsif n = 1
         then declare
    	        min,max : integer := 0;
	        d : Link_to_Vector := Head_Of(al(al'first));
              begin
	        Min_Max(al(al'first),d'first,min,max);
	        return (max - min);
              end;
         elsif All_Equal(al)
	     then return Volume(n,al(al'first));
	     else return Minkowski_Sum(n,al);
    end if;
  end Mixed_Volume;

-- MIXED VOLUME COMPUTATIONS WITH TREE OF USEFUL DIRECTIONS :

  function Two_Terms_Mixed_Volume ( n : natural; al : Array_of_Lists;
 				    tv : Tree_of_Vectors )
		       	          return natural is

  -- DESCRIPTION :
  --   returns the mixed volume of the polytopes generated by the
  --   points in al, where the first polytope is generated by only
  --   two points.

    first,second : Link_to_Vector;
    l : List renames al(al'first);
    res : natural;

  begin
    first := Head_Of(l);
    second := Head_Of(Tail_Of(l));
    declare
      d : Vector(first'range);
      piv : integer := 0;
    begin
      for i in d'range loop
	d(i) := first(i) - second(i);
	if (piv = 0) and then (d(i) /= 0)
	 then piv := i;
        end if;
      end loop;
      if piv = 0
       then return 0;
       else if d(piv) < 0
	     then Min_Vector(d);
            end if;
	    declare
	      t : Transfo := Rotate(d,piv);
	      tr_al : Array_of_Lists(al'first..(al'last-1));
              degen : boolean := false;
            begin
              Apply(t,d);
	      for i in tr_al'range loop
		tr_al(i) := Transform_and_Reduce(t,piv,al(i+1));
                Remove_Duplicates(tr_al(i));
                degen := (Length_Of(tr_al(i)) <= 1);
                exit when degen;
              end loop;
              Clear(t);
              if not degen
	       then res := d(piv)*Mixed_Volume(n-1,tr_al,tv);
               else res := 0;
              end if;
              Deep_Clear(tr_al);
            end;
      end if;
    end;
    return res;
  end Two_Terms_Mixed_Volume;
     
  function Minkowski_Sum ( n : natural; al : Array_of_Lists;
                           tv : Tree_of_Vectors ) return natural is

  -- DESCRIPTION :
  --   Computes the mixed volume of the polytope generated
  --   by the points in al, where n > 1.
  --   The tree of degrees is not empty.

    res,m : natural;
    al1 : List;
    wal : Array_of_Lists(al'range) := Interchange2(al);
    tmp : Tree_of_Vectors;
    perm,mix : Link_to_Vector;

  begin
    m := Length_Of(wal(wal'first));
    if m = 2
     then return Two_Terms_Mixed_Volume(n,wal,tv);
     elsif m > 2
         then
           Mixture(al,perm,mix);
           wal := Permute(perm.all,al);
           declare
	     shiftvec : Link_to_Vector;
             nulvec : Vector(1..n) := (1..n => 0);
	     shifted : boolean;
           begin
            -- SHIFT OF THE FIRST LIST ( then all pv >= 0) :
             if not Is_In(wal(wal'first),nulvec)
              then shiftvec := Graded_Max(wal(wal'first));
      	           al1 := Shift(wal(wal'first),shiftvec); shifted := true;
              else al1 := wal(wal'first);                 shifted := false;
             end if;
             Move_to_Front(al1,nulvec);
             wal(wal'first) := al1;
	    -- COMPUTING THE MIXED VOLUME :
             tmp := tv;  res := 0;
	     while not Is_Null(tmp) loop
	       declare
	         nd : node := Head_Of(tmp);
	         v : Link_to_Vector := nd.d;
	         i : integer := Pivot(v);
	         pv : integer := Maximal_Support(al1,v.all);
	         pal : Array_of_Lists(al'first..(al'last-1));
	         degen : boolean;
               begin
	         Projection(wal,v.all,i,pal,degen);
	         if (nd.ltv = null) or else Is_Null(nd.ltv.all)
                  then res := res + pv*Mixed_Volume(n-1,pal);
                  else res := res + pv*Mixed_Volume(n-1,pal,nd.ltv.all);
                 end if;
                 Deep_Clear(pal);
               end;
               tmp := Tail_Of(tmp);
             end loop;
            -- CLEANING UP :
             Clear(perm); Clear(mix);
             if shifted 
	      then Clear(shiftvec); Deep_Clear(al1);
             end if;
             return res;
           end;
         else -- m < 2
           return 0;
    end if;
  end Minkowski_Sum;

  function Mixed_Volume ( n : natural; al : Array_of_Lists;
			  tv : Tree_of_Vectors ) return natural is
  begin
    if (n = 0) or Is_Null(al(al'first))
     then return 0;
     elsif n = 1
         then declare
    	        min,max : integer := 0;
	        d : Link_to_Vector := Head_Of(al(al'first));
              begin
	        Min_Max(al(al'first),d'first,min,max);
	        return (max - min);
              end;
         elsif All_Equal(al)
	     then return Volume(n,al(al'first),tv);
	     elsif not Is_Null(tv)
		 then return Minkowski_Sum(n,al,tv);
		 else return Minkowski_Sum(n,al);
    end if;
  end Mixed_Volume;

-- MIXED VOLUME COMPUTATIONS WITH CREATION OF TREE OF USEFUL DIRECTIONS :

  procedure Two_Terms_Mixed_Vol 
             ( n : in natural; al : in Array_of_Lists;
	       tv : in out Tree_of_Vectors; mv : out natural ) is

  -- DESCRIPTION :
  --   returns the mixed volume of the polytopes generated by the
  --   points in al, where the first polytope is generated by only
  --   two points.

    first,second : Link_to_Vector;
    l : List renames al(al'first);

  begin
    first := Head_Of(l);
    second := Head_Of(Tail_Of(l));
    declare
      d : Vector(first'range);
      piv : integer := 0;
    begin
      for i in d'range loop
	d(i) := first(i) - second(i);
	if (piv = 0) and then (d(i) /= 0)
	 then piv := i;
        end if;
      end loop;
      if piv = 0
       then mv := 0;
       else if d(piv) < 0
	     then Min_Vector(d);
            end if;
	    declare
	      t : Transfo := Rotate(d,piv);
	      tr_al : Array_of_Lists(al'first..(al'last-1));
	      mvl : natural;
              degen : boolean := false;
            begin
              Apply(t,d);
	      for i in tr_al'range loop
		tr_al(i) := Transform_and_Reduce(t,piv,al(i+1));
                Remove_Duplicates(tr_al(i));
                degen := (Length_Of(tr_al(i)) <= 1);
                exit when degen;
              end loop;
              Clear(t);
              if not degen
	       then Mixed_Volume(n-1,tr_al,tv,mvl); mv := d(piv)*mvl;
               else mv := 0;
              end if;
              Deep_Clear(tr_al);
            end;
      end if;
    end;
  end Two_Terms_Mixed_Vol;
     
  procedure Minkowski_Sum ( n : in natural; al : in Array_of_Lists;
                            tv : in out Tree_of_Vectors; mv : out natural ) is

  -- DESCRIPTION :
  --   Computes the mixed volume of the polytope generated
  --   by the points in al, where n > 1.

    res,m,ptslen : natural;
    al1 : List;
    typ,ind : Vector(1..n);
    perm,mix : Link_to_Vector;
    wal : Array_of_Lists(al'range) := Interchange2(al);

    procedure Update ( v : in Vector; i : in integer;
		       added : in out boolean ) is

    -- DESCRIPTION :
    --   This procedure computes the support of the first list
    --   in the direction v; if this support is not zero,
    --   the projection will be computed.
    --   The parameter added becomes true if v has been added to vl.

      pal : Array_of_Lists(al'first..(al'last-1));
      pv : integer;
      degen : boolean;

    begin
      if not Is_In(tv,v)
       then pv := Maximal_Support(al1,v);
            if pv > 0
             then Projection(wal,v,i,pal,degen);
                  if not degen
                   then declare
                          tmp : Tree_of_Vectors;
                          mv : natural;
                        begin
                          Mixed_Volume(n-1,pal,tmp,mv);
                          if mv = 0
                           then Clear(tmp);
                           else res := res + pv*mv;
                                declare
                                  nd : node;
                                begin
                                  nd.d := new Integer_Vectors.Vector'(v);
                                  if not Is_Null(tmp)
                                   then nd.ltv := new Tree_of_Vectors'(tmp);
                                  end if;
                                  Construct(nd,tv);
                                end;
                                added := true;
                          end if;
                        end;
                        Deep_Clear(pal);
                  end if;
            end if;
      end if;
    end Update;

    procedure Enumerate_Facets ( lpts : in Array_of_Lists; len : in natural ) is

      pts : Integer_Vectors_of_Vectors.Vector(1..len);
      cnt : integer;

      procedure Compute_Facet
                   ( facet : in Integer_Vectors_of_Vectors.Vector;
                     cont : out boolean ) is

        v,w : Vector(1..n);
        i : integer;
        added : boolean;
      begin
        v := Facet_Normal(n,facet,pts);
        i := Pivot(v);
        if i <= v'last
         then added := false;
              Update(v,i,added);
              if not added
               then Min_Vector(v); w := v;
               else w := -v; added := false;
              end if;
              Update(w,i,added);
        end if;
        cont := true;
      end Compute_Facet;

      procedure Compute_Facets is
        new Face_Enumerator_of_Sum.Enumerate_Faces_of_Sum(Compute_Facet);
       -- new Face_Enumerators.Enumerate_Faces_of_Sum(Compute_Facet);

    begin
      pts(ind(1)..ind(2)-1) := Shallow_Create(Tail_Of(al1));
      cnt := lpts'first + mix(mix'first);
      for i in mix'first+1..mix'last loop
        if i < mix'last
         then pts(ind(i)..ind(i+1)-1) := Shallow_Create(lpts(cnt));
              cnt := cnt + mix(i);
         else pts(ind(i)..len) := Shallow_Create(lpts(cnt));
        end if;
      end loop;
      Compute_Facets(ind(mix'range),typ(mix'range),n-1,pts);
    end Enumerate_Facets;

  begin
    m := Length_Of(wal(wal'first));
    if m = 2
     then Two_Terms_Mixed_Vol(n,wal,tv,mv);
     elsif m > 2
         then
           Mixture(al,perm,mix);
           wal := Permute(perm.all,al);
           declare
	     shiftvec : Link_to_Vector;
             nulvec : Vector(1..n) := (1..n => 0);
             shifted : boolean;
             cnt : integer;
           begin
            -- SHIFT OF THE FIRST LIST ( then all pv >= 0) :
             if not Is_In(wal(wal'first),nulvec)
              then shiftvec := Graded_Max(wal(wal'first));
                   al1 := Shift(wal(wal'first),shiftvec); shifted := true;
              else al1 := wal(wal'first);                 shifted := false;
             end if;
            -- ENUMERATE FACES OF SUM OF THE RIGHT TYPE :
             Move_to_Front(al1,nulvec);
             wal(wal'first) := al1;
             res := 0;
             typ(1) := mix(mix'first)-1;
             typ(2) := mix(mix'first+1); ind(1) := 1;
             ind(2) := ind(1) + Length_Of(al1) - 1;  -- skip null vector
             cnt := wal'first + mix(mix'first);
             for i in mix'first+2..mix'last loop
               typ(i) := mix(i);
               ind(i) := ind(i-1) + Length_Of(wal(cnt));
               cnt := cnt + mix(i-1);
             end loop;
             ptslen := ind(mix'last) + Length_Of(wal(cnt)) - 1;
             Enumerate_Facets(wal,ptslen);
            -- CLEANING UP :
             Clear(perm); Clear(mix);
             if shifted 
              then Clear(shiftvec); Deep_Clear(al1);
             end if;
             mv := res;
           end;
         else -- m < 2
           mv := 0;
    end if;
  end Minkowski_Sum;

  procedure Mixed_Volume ( n : in natural; al : in Array_of_Lists;
			   tv : in out Tree_of_Vectors; mv : out natural ) is
  begin
    if (n = 0) or Is_Null(al(al'first))
     then mv := 0;
     elsif n = 1
         then declare
    	        min,max : integer := 0;
	        d : Link_to_Vector := Head_Of(al(al'first));
              begin
	        Min_Max(al(al'first),d'first,min,max);
	        mv := max - min;
              end;
         elsif All_Equal(al)
	     then Volume(n,al(al'first),tv,mv);
	     else Minkowski_Sum(n,al,tv,mv);
    end if;
  end Mixed_Volume;

end Volumes;