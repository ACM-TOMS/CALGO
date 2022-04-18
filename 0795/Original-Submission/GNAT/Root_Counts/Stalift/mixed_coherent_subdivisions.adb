with Integer_Faces_of_Polytope;         use Integer_Faces_of_Polytope;
with Integer_Lifting_Functions;         use Integer_Lifting_Functions;
with Integer_Pruning_Methods;           use Integer_Pruning_Methods;

package body Mixed_Coherent_Subdivisions is

-- a polynomial system as lifting function :

  function Mixed_Coherent_Subdivision
               ( n : natural; mix : Vector; points : Array_of_Lists; 
                 lift : Poly_Sys ) return Mixed_Subdivision is

    res : Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);
    nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);

  begin
    Mixed_Coherent_Subdivision(n,mix,points,lift,lifted,nbsucc,nbfail,res);
    Deep_Clear(lifted);
    return res;
  end Mixed_Coherent_Subdivision;

  procedure Mixed_Coherent_Subdivision
               ( n : in natural; mix : in Vector; points : in Array_of_Lists;
                 lift : in Poly_Sys; lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision ) is

    fa : Array_of_Faces(mix'range);
    index : natural := points'first;

  begin
    for k in lifted'range loop                        -- generate lower faces
      lifted(k) := Polynomial_Lift(lift(k),points(index));
      fa(k) := Create_Lower(mix(k),n+1,lifted(k));
      index := index + mix(k);
    end loop;
    Create_CS(n,mix,fa,lifted,nbsucc,nbfail,mixsub);  -- prune for mixed cells
    Shallow_Clear(fa);
  end Mixed_Coherent_Subdivision;

-- a user-defined lifting function :

  function Mixed_Coherent_Subdivision
               ( n : natural; mix : Vector; points : Array_of_Lists;
                 linear : boolean; lift : Integer_Vectors_of_Vectors.Vector )
               return Mixed_Subdivision is

    res : Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);
    nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);

  begin
    Mixed_Coherent_Subdivision
      (n,mix,points,linear,lift,lifted,nbsucc,nbfail,res);
    Deep_Clear(lifted);
    return res;
  end Mixed_Coherent_Subdivision;

  procedure Mixed_Coherent_Subdivision
               ( n : in natural; mix : in Vector; points : in Array_of_Lists;
                 linear : in boolean;
                 lift : in Integer_Vectors_of_Vectors.Vector;
                 lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision ) is

    fa : Array_of_Faces(mix'range);
    index : natural := points'first;

  begin
    for k in lifted'range loop                          -- compute lower faces
      if linear
       then lifted(k) := Linear_Lift(lift(k).all,points(index));
       else lifted(k) := Point_Lift(lift(k).all,points(index));
      end if;
      fa(k) := Create_Lower(mix(k),n+1,lifted(k));
      index := index + mix(k);
    end loop;
    Create_CS(n,mix,fa,lifted,nbsucc,nbfail,mixsub);  -- prune for mixed cells
    Shallow_Clear(fa);
  end Mixed_Coherent_Subdivision;

-- a randomly generated lifting function :

  function Mixed_Coherent_Subdivision
               ( n : natural; mix : Vector; points : Array_of_Lists;
                 linear : boolean; low,upp : Vector )
               return Mixed_Subdivision is

    res : Mixed_Subdivision;
    lifted : Array_of_Lists(mix'range);
    nbsucc,nbfail : Float_Vectors.Vector(mix'range) := (mix'range => 0.0);

  begin
    Mixed_Coherent_Subdivision
      (n,mix,points,linear,low,upp,lifted,nbsucc,nbfail,res);
    Deep_Clear(lifted);
    return res;
  end Mixed_Coherent_Subdivision;

  procedure Mixed_Coherent_Subdivision
               ( n : in natural; mix : in Vector; points : in Array_of_Lists;
                 linear : in boolean; low,upp : in Vector;
                 lifted : in out Array_of_Lists;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision ) is

    fa : Array_of_Faces(mix'range);
    index : natural := points'first;

  begin
    for k in lifted'range loop                          -- compute lower faces
      if linear
       then lifted(k) := Random_Linear_Lift(low(k),upp(k),points(index));
       else lifted(k) := Random_Lift(low(k),upp(k),points(index));
      end if;
      fa(k) := Create_Lower(mix(k),n+1,lifted(k));
      index := index + mix(k);
    end loop;
    Create_CS(n,mix,fa,lifted,nbsucc,nbfail,mixsub);  -- prune for mixed cells
    Shallow_Clear(fa);
  end Mixed_Coherent_Subdivision;

end Mixed_Coherent_Subdivisions;
