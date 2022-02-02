with Integer_Vectors_of_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Dynamic_Triangulations;             use Dynamic_Triangulations;
with Unfolding_Subdivisions;             use Unfolding_Subdivisions;

package body Triangulations_and_Subdivisions is

-- REFINEMENT ROUTINES :

  procedure Refine ( n : in natural; mic : in out Mixed_Cell ) is

  -- NOTE :
  --   Dynamic lifting will be applied with standard settings,
  --   under the assumption that there are only few points in the cell.

    support : List := Reduce(mic.pts(1),n+1);
    t : Triangulation;
    lifted,lifted_last : List;

  begin
    Dynamic_Lifting(support,false,true,0,lifted,lifted_last,t);
    mic.sub := new Mixed_Subdivision'(Deep_Create(n,t));
    Deep_Clear(lifted); Clear(t); 
      -- pity that Shallow_Clear(t) is not yet possible ...
  end Refine;

  procedure Refine ( n : in natural; mixsub : in out Mixed_Subdivision ) is

  -- NOTE :
  --   Refines the mixed subdivision, under the safe assumption that 
  --   there is only one support set to deal with.

    res,res_last : Mixed_Subdivision;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      if Length_Of(mic.pts(1)) > n+1
       then Refine(n,mic);
      end if;
      Append(res,res_last,mic);
      tmp := Tail_Of(tmp);
    end loop;
    mixsub := res;
  end Refine;

-- TARGET PROCEDURES :

  function Deep_Create ( n : natural; s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;
    ver : constant Integer_Vectors_of_Vectors.Vector := Vertices(s);

  begin
    res.nor := new Integer_Vectors.Vector'(Normal(s));
    res.pts := new Array_of_Lists(1..1);
    res.pts(1) := Deep_Create(ver);
    return res;
  end Deep_Create;

  function Shallow_Create ( n : natural; s : Simplex ) return Mixed_Cell is

    res : Mixed_Cell;
    ver : constant Integer_Vectors_of_Vectors.Vector := Vertices(s);

  begin
    res.nor := new Integer_Vectors.Vector'(Normal(s));
    res.pts := new Array_of_Lists(1..1);
    res.pts(1) := Shallow_Create(ver);
    return res;
  end Shallow_Create;

  function Deep_Create ( n : natural; t : Triangulation )
                       return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Deep_Create(n,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Deep_Create;

  function Shallow_Create ( n : natural; t : Triangulation )
                          return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;

  begin
    while not Is_Null(tmp) loop
      Append(res,res_last,Shallow_Create(n,Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Shallow_Create;

  function Deep_Create ( n : natural; flatnor : Vector; t : Triangulation )
                       return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;
    s : Simplex;

  begin
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      exit when (flatnor = Normal(s));
      Append(res,res_last,Deep_Create(n,s));
      tmp := Tail_Of(tmp);
    end loop;
    res := Merge(res);             -- merge cells with same inner normal
    Refine(n,res);                 -- refine the non-fine cells
    return res;
  end Deep_Create;

  function Shallow_Create ( n : natural; flatnor : Vector; t : Triangulation )
                          return Mixed_Subdivision is

    res,res_last : Mixed_Subdivision;
    tmp : Triangulation := t;
    s : Simplex;

  begin
    while not Is_Null(tmp) loop
      s := Head_Of(tmp);
      exit when (flatnor = Normal(s));
      Append(res,res_last,Shallow_Create(n,s));
      tmp := Tail_Of(tmp);
    end loop;
    res := Merge(res);             -- merge cells with same inner normal
    Refine(n,res);                 -- refine the non-fine cells
    return res;
  end Shallow_Create;

  function Non_Flat_Deep_Create ( n : natural; t : Triangulation )
                                return Mixed_Subdivision is

    flatnor : Vector(1..n+1) := (1..n+1 => 0);

  begin
    flatnor(n+1) := 1;
    return Deep_Create(n,flatnor,t);
  end Non_Flat_Deep_Create;

  function Non_Flat_Shallow_Create ( n : natural; t : Triangulation )
                                   return Mixed_Subdivision is

    flatnor : Vector(1..n+1) := (1..n+1 => 0);

  begin
    flatnor(n+1) := 1;
    return Shallow_Create(n,flatnor,t);
  end Non_Flat_Shallow_Create;

  function Deep_Create ( n : natural; mixsub : Mixed_Subdivision )
                       return Triangulation is

    res : Triangulation;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      declare
        v : Integer_Vectors_of_Vectors.Vector(0..n);
        tmppts : List := mic.pts(mic.pts'first);
        s : Simplex;
      begin
        for i in v'range loop
          v(i) := new Integer_Vectors.Vector'(Head_Of(tmppts).all);
          tmppts := Tail_Of(tmppts);
          exit when Is_Null(tmppts);
        end loop;
        s := Create(v);
        Construct(s,res);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Connect(res);
    return res;
  end Deep_Create;

  function Shallow_Create ( n : natural; mixsub : Mixed_Subdivision )
                          return Triangulation is

    res : Triangulation;
    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      declare
        v : Integer_Vectors_of_Vectors.Vector(0..n);
        tmppts : List := mic.pts(mic.pts'first);
        s : Simplex;
      begin
        for i in v'range loop
          v(i) := Head_Of(tmppts);
          tmppts := Tail_Of(tmppts);
          exit when Is_Null(tmppts);
        end loop;
        s := Create(v);
        Construct(s,res);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Connect(res);
    return res;
  end Shallow_Create;

end Triangulations_and_Subdivisions;
