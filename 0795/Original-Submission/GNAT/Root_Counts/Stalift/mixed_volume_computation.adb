with integer_io,Integer_Vectors_io;      use integer_io,Integer_Vectors_io;
with Arrays_of_Integer_Vector_Lists_io;  use Arrays_of_Integer_Vector_Lists_io;

with Float_Vectors;
with Integer_Matrices;                   use Integer_Matrices;
with Integer_Linear_System_Solvers;      use Integer_Linear_System_Solvers;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Mixed_Coherent_Subdivisions;        use Mixed_Coherent_Subdivisions;

package body Mixed_Volume_Computation is

-- AUXILIAIRY OUTPUT ROUTINES :

  procedure put ( file : in file_type; points : in Array_of_Lists;
                  n : in natural; mix : in Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural ) is
  begin
    new_line(file);
    put_line(file,"THE LIFTED SUPPORTS :");
    new_line(file);
    put(file,points);
    new_line(file);
    put_line(file,"THE MIXED SUBDIVISION :");
    new_line(file);
    put(file,n,mix,mixsub,mv);
  end put;

  procedure Sort ( supports : in out Array_of_Lists; k,nb,n : in natural;
                   mxt,perm : in out Vector ) is

  -- DESCRIPTION :
  --   Auxiliary operation for Compute_Mixture.
  --   Compares the kth support with the following supports.
  --   Already nb different supports have been found.

  begin
    for l in (k+1)..n loop
      if Is_Equal(Supports(k),Supports(l))
       then if l /= k + mxt(nb)
             then declare
                    pos : natural := k + mxt(nb);
                    tmpdl : List := supports(l);
                    tmppos : natural;
                  begin
                    supports(l) := supports(pos);
                    supports(pos) := tmpdl;
                    tmppos := perm(l);
                    perm(l) := perm(pos);
                    perm(pos) := tmppos;
                  end;
            end if;
            mxt(nb) := mxt(nb) + 1;
      end if;
    end loop;
  end Sort;

-- TARGET ROUTINES :

  procedure Compute_Mixture ( supports : in out Array_of_Lists;
			      mix,perms : out Link_to_Vector ) is

    n : constant natural := supports'last;
    cnt : natural := 0;              -- counts the number of different supports
    mxt : Vector(supports'range)     -- counts the number of occurrencies
        := (supports'range => 1);
    perm : Link_to_Vector            -- keeps track of the permutations
         := new Integer_Vectors.Vector(supports'range);
    index : natural := supports'first;

  begin
    for k in perm'range loop
      perm(k) := k;
    end loop;
    while index <= supports'last loop
      cnt := cnt + 1;
      Sort(supports,index,cnt,n,mxt,perm.all);
      index := index + mxt(cnt);
    end loop;
    mix := new Integer_Vectors.Vector'(mxt(mxt'first..cnt));
    perms := perm;
  end Compute_Mixture;

  function Compute_Index ( k : natural; mix : Vector ) return natural is

  -- DESCRIPTION :
  --   Returns the index of k w.r.t. to the type of mixture.

    index : natural := mix(mix'first);

  begin
    if k <= index
     then return mix'first;
     else for l in (mix'first+1)..mix'last loop
            index := index + mix(l);
            if k <= index
             then return l;
            end if;
          end loop;
          return mix'last;
    end if;
  end Compute_Index;

  function Compute_Permutation ( n : natural; mix : Vector;
                                 supports : Array_of_Lists )
                               return Link_to_Vector is

    perms : Link_to_Vector := new Vector(1..n);

  begin
    for k in perms'range loop
      perms(k) := k;
    end loop;
    return perms;
  end Compute_Permutation;

  function Permute ( p : Poly_Sys; perm : Link_to_Vector ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for k in p'range loop
      res(k) := p(perm(k));
    end loop;
    return res;
  end Permute;

  function Permute ( supports : Array_of_Lists; perm : Link_to_Vector )
                   return Array_of_Lists is

    res : Array_of_Lists(supports'range);

  begin
    for k in supports'range loop
      res(k) := supports(perm(k));
    end loop;
    return res;
  end Permute;

  function Typed_Lists ( mix : Vector; points : Array_of_Lists )
                       return Array_of_Lists is

    res : Array_of_Lists(mix'range);
    ind : natural := res'first;

  begin
    for i in mix'range loop
      res(i) := points(ind);
      ind := ind + mix(i);
    end loop;
    return res;
  end Typed_Lists;

-- MIXED VOLUME COMPUTATIONS BASED ON SUBDIVISIONS :

-- AUXILIARIES :

  function Is_Fine ( mix : Vector; mic : Mixed_Cell ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the mixed volume can be computed by a determinant.

    fine : boolean := true;

  begin
    for k in mic.pts'range loop
      fine := (Length_Of(mic.pts(k)) = mix(k) + 1);
      exit when not fine;
    end loop;
    return fine;
  end Is_Fine;

  function Reduced_Supports ( n : natural; mix : Vector; mic : Mixed_Cell )
                            return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the supports of the cell without the lifting values.

    res : Array_of_Lists(1..n);
    cnt : natural := 1;

  begin
    for k in mic.pts'range loop
      res(cnt) := Reduce(mic.pts(k),n+1);
      for l in 1..mix(k)-1 loop
        Copy(res(cnt),res(cnt+l));
      end loop;
      cnt := cnt + mix(k);
    end loop;
    return res;
  end Reduced_Supports;

  function Fine_Mixed_Volume ( n : natural; mix : Vector; mic : Mixed_Cell )
                             return natural is

  -- DESCRIPTION :
  --   Computes the mixed volume for a cell that is fine mixed.

  -- REQUIRED : Fine(mix,mic).

    res,count : natural;
    mat : matrix(1..n,1..n);
    detmat : integer;
    tmp : List;
    sh,pt : Link_to_Vector;

  begin
    count := 1;
    for k in mic.pts'range loop
      sh := Head_Of(mic.pts(k));
      tmp := Tail_Of(mic.pts(k));
      while not Is_Null(tmp) loop
        pt := Head_Of(tmp);
        for j in 1..n loop
          mat(count,j) := pt(j) - sh(j);
        end loop;
        tmp := Tail_Of(tmp);
        count := count + 1;
      end loop;
    end loop;
    detmat := Det(mat);
    if detmat >= 0
     then res := detmat;
     else res := -detmat;
    end if;
    return res;
  end Fine_Mixed_Volume;

  function Mixed_Volume ( n : natural; mix : Vector;
			  mic : Mixed_Cell ) return natural is

  -- ALGORITHM :
  --   First check if the cell has a refinement, if so, then use it,
  --   if not, then check if the cell is fine mixed.
  --   If the cell is fine mixed, only a determinant needs to be computed,
  --   otherwise the cell will be refined.
   
    res : natural;

  begin
    if (mic.sub /= null) and then not Is_Null(mic.sub.all)
     then res := Mixed_Volume_Computation.Mixed_Volume(n,mix,mic.sub.all);
     elsif Is_Fine(mix,mic)
         then res := Fine_Mixed_Volume(n,mix,mic);
         else declare
                rcell : Array_of_Lists(1..n) := Reduced_Supports(n,mix,mic);
              begin
                res := Mixed_Volume_Computation.Mixed_Volume(n,rcell);
                Deep_Clear(rcell);
              end;
    end if;
    return res;
  end Mixed_Volume;

  function Mixed_Volume ( n : natural; mix : Vector;
			  mixsub : Mixed_Subdivision ) return natural is

    res : natural := 0;
    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      res := res + Mixed_Volume(n,mix,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Mixed_Volume;

  procedure Mixed_Volume ( n : in natural; mix : in Vector;
                           mic : in out Mixed_Cell; mv : out natural ) is
  begin
    if (mic.sub /= null) and then not Is_Null(mic.sub.all)
     then mv := Mixed_Volume_Computation.Mixed_Volume(n,mix,mic.sub.all);
     elsif Is_Fine(mix,mic)
         then mv := Fine_Mixed_Volume(n,mix,mic);
         else -- NOTE : keep the same type of mixture!
           declare
             rcell : Array_of_Lists(1..n) := Reduced_Supports(n,mix,mic);
             lifted : Array_of_Lists(mix'range);
             mixsub : Mixed_Subdivision;
           begin
             Mixed_Volume_Computation.Mixed_Volume
               (n,mix,rcell,lifted,mixsub,mv);
             mic.sub := new Mixed_Subdivision'(mixsub);
             Deep_Clear(rcell);  Deep_Clear(lifted);
           end;
    end if;
  end Mixed_Volume;

  procedure Mixed_Volume ( n : in natural; mix : in Vector;
                           mixsub : in out Mixed_Subdivision;
                           mv : out natural ) is

    res : natural := 0;
    tmp : Mixed_Subdivision := mixsub;

  begin
    while not Is_Null(tmp) loop
      declare
        mic : Mixed_Cell := Head_Of(tmp);
        mmv : natural;
      begin
        Mixed_Volume(n,mix,mic,mmv);
        Set_Head(tmp,mic);
        res := res + mmv;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    mv := res;
  end Mixed_Volume;

-- MIXED VOLUME COMPUTATIONS BASED ON SUPPORTS :

  function Mixed_Volume ( n : natural; supports : Array_of_Lists )
			return natural is
    mv : natural;
    mix,perm : Link_to_Vector;
    permsupp : Array_of_Lists(supports'range);

  begin
    Copy(supports,permsupp);
    Compute_Mixture(permsupp,mix,perm);
    mv := Mixed_Volume(n,mix.all,permsupp);
    Clear(mix); Clear(perm); Deep_Clear(permsupp);
    return mv;
  end Mixed_Volume;

  function Mixed_Volume ( file : file_type; n : natural;
                          supports : Array_of_Lists ) return natural is
    mv : natural;
    mix,perm : Link_to_Vector;
    permsupp : Array_of_Lists(supports'range);

  begin
    Copy(supports,permsupp);
    Compute_Mixture(permsupp,mix,perm);
    mv := Mixed_Volume(file,n,mix.all,permsupp);
    Clear(mix); Clear(perm); Deep_Clear(permsupp);
    return mv;
  end Mixed_Volume;

  function Mixed_Volume ( n : natural; mix : Vector;
                          supports : Array_of_Lists ) return natural is

    res : natural;
    mixsub : Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);

  begin
    Mixed_Volume_Computation.Mixed_Volume(n,mix,supports,lifted,mixsub,res);
    Deep_Clear(lifted);
    Shallow_Clear(mixsub);
    return res;
  end Mixed_Volume;

  procedure Mixed_Volume 
              ( n : in natural; mix : in Vector; supports : in Array_of_Lists;
                lifted : out Array_of_Lists;
                mixsub : out Mixed_Subdivision; mv : out natural ) is

    low : constant Vector := (mix'range => 0);
    upp : constant Vector := Adaptive_Lifting(supports);
    nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);
    liftsupp : Array_of_Lists(mix'range);
    sub : Mixed_Subdivision;

  begin
    Mixed_Coherent_Subdivision(n,mix,supports,false,low,upp,liftsupp,
                               nbsucc,nbfail,sub);
    Mixed_Volume(n,mix,sub,mv);
    lifted := liftsupp;
    mixsub := sub;
  end Mixed_Volume;

  function Mixed_Volume ( file : file_type; n : natural; mix : Vector;
                          supports : Array_of_Lists ) return natural is
  
    res : natural;
    mixsub : Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);

  begin
    Mixed_Volume_Computation.Mixed_Volume
      (file,n,mix,supports,lifted,mixsub,res);
    Deep_Clear(lifted);
    Shallow_Clear(mixsub);
    return res;
  end Mixed_Volume;

  procedure Mixed_Volume
              ( file : in file_type; n : in natural;
                mix : in Vector; supports : in Array_of_Lists;
                lifted : out Array_of_Lists;
                mixsub : out Mixed_Subdivision; mv : out natural ) is

    low : constant Vector := (mix'range => 0);
    upp : constant Vector := Adaptive_Lifting(supports);
    sub : Mixed_Subdivision;
    nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);
    liftsupp : Array_of_Lists(mix'range);

  begin
    Mixed_Coherent_Subdivision
      (n,mix,supports,false,low,upp,liftsupp,nbsucc,nbfail,sub);
    lifted := liftsupp;
    mixsub := sub;
    put(file,liftsupp,n,mix,sub,mv);
  end Mixed_Volume;

end Mixed_Volume_Computation;
