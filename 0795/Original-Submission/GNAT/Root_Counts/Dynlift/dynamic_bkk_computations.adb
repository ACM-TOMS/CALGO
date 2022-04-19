with Integer_Vectors,Float_Vectors;      use Integer_Vectors;
with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Power_Lists,Vertices;               use Power_Lists,Vertices;

with Triangulations,Triangulations_io;   use Triangulations,Triangulations_io;
with Dynamic_Triangulations;             use Dynamic_Triangulations;
with Integer_Faces_of_Polytope;          use Integer_Faces_of_Polytope;
with Integer_Mixed_Subdivisions;         use Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;      use Integer_Mixed_Subdivisions_io;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Dynamic_Mixed_Subdivisions;         use Dynamic_Mixed_Subdivisions;
with Dynamic_Polyhedral_Continuation;    use Dynamic_Polyhedral_Continuation;

package body Dynamic_BKK_Bound_Computations is

  function BKK_by_Dynamic_Lifting ( p : Poly_Sys ) return natural is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    mix,perms : Link_to_Vector;
    res : natural;

  begin
    Compute_Mixture(supports,mix,perms);
    if mix'first = mix'last
     then declare
            support : List := supports(supports'first);
            vert : List := Vertex_Points(support);
            lifted,lifted_last : List;
            t : Triangulation;
          begin
            Dynamic_Lifting(vert,false,false,0,lifted,lifted_last,t);
            res := Volume(t);
            Deep_Clear(vert); Deep_Clear(lifted); Clear(t);
          end;
     else declare
            verpts : Array_of_Lists(mix'range);
            cnt : natural := verpts'first;
            mixsub : Mixed_Subdivision;
            fs : Face_Structures(mix'range);
            nbsucc,nbfail : Float_Vectors.Vector(mix'range)
                          := (mix'range => 0.0);
          begin
            for i in mix'range loop
              verpts(cnt) := Vertex_Points(supports(cnt));
              cnt := cnt + mix(i);
            end loop;
            Dynamic_Lifting(n,mix.all,verpts,false,false,false,0,mixsub,fs,
                            nbsucc,nbfail);
            res := Mixed_Volume(n,mix.all,mixsub);
            Deep_Clear(verpts);
            for i in fs'range loop
              Clear(fs(i).t); Deep_Clear(fs(i).l); Deep_Clear(fs(i).f);
            end loop;
          end;
    end if;
    return res;
  end BKK_by_Dynamic_Lifting;

  function BKK_by_Dynamic_Lifting ( file : file_type; p : Poly_Sys )
                                  return natural is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    mix,perms : Link_to_Vector;
    res : natural;

  begin
    Compute_Mixture(supports,mix,perms);
    if mix'first = mix'last
     then declare
            support : List := supports(supports'first);
            vert : List := Vertex_Points(support);
            lifted,lifted_last : List;
            t : Triangulation;
          begin
            Dynamic_Lifting(vert,false,false,0,lifted,lifted_last,t);
            put(file,n,t,res);
            Deep_Clear(vert); Deep_Clear(lifted); Clear(t);
          end;
     else declare
            verpts : Array_of_Lists(mix'range);
            cnt : natural := verpts'first;
            mixsub : Mixed_Subdivision;
            fs : Face_Structures(mix'range);
            nbsucc,nbfail : Float_Vectors.Vector(mix'range)
                          := (mix'range => 0.0);
          begin
            for i in mix'range loop
              verpts(cnt) := Vertex_Points(supports(cnt));
              cnt := cnt + mix(i);
            end loop;
            Dynamic_Lifting(n,mix.all,verpts,false,false,false,0,mixsub,fs,
                            nbsucc,nbfail);
            put(file,n,mix.all,mixsub,res);
            Deep_Clear(verpts);
            for i in fs'range loop
              Clear(fs(i).t); Deep_Clear(fs(i).l); Deep_Clear(fs(i).f);
            end loop;
          end;
    end if;
    return res;
  end BKK_by_Dynamic_Lifting;

  function Solve_by_Dynamic_Lifting ( p : Poly_Sys ) return Solution_List is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    mix,perms : Link_to_Vector;
    sols : Solution_List;
    file : file_type;

  begin
    Create(file,out_file,"/tmp/brol");
    Compute_Mixture(supports,mix,perms);
    if mix'first = mix'last
     then declare
            support : List := supports(supports'first);
            vert : List := Vertex_Points(support);
            lifted,lifted_last : List;
            t : Triangulation;
          begin
            Dynamic_Unmixed_Solve(file,n,vert,false,false,0,lifted,lifted_last,
                                  t,p,sols);
           -- Deep_Clear(vert); Deep_Clear(lifted); Clear(t);
          end;
     else declare
            verpts : Array_of_Lists(mix'range);
            cnt : natural := verpts'first;
            mixsub : Mixed_Subdivision;
            fs : Face_Structures(mix'range);
            nbsucc,nbfail : Float_Vectors.Vector(mix'range)
                          := (mix'range => 0.0);
          begin
            for i in mix'range loop
              verpts(cnt) := Vertex_Points(supports(cnt));
              cnt := cnt + mix(i);
            end loop;
            Dynamic_Mixed_Solve(file,n,mix.all,supports,false,false,false,0,
                                mixsub,fs,nbsucc,nbfail,p,sols);
           -- Deep_Clear(verpts);
           -- for i in fs'range loop
           --   Clear(fs(i).t); Deep_Clear(fs(i).l); Deep_Clear(fs(i).f);
           -- end loop;
          end; 
    end if;
    Close(file);
    return sols;
  end Solve_by_Dynamic_Lifting;

  function Solve_by_Dynamic_Lifting ( file : file_type; p : Poly_Sys )
                                    return Solution_List is
    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    mix,perms : Link_to_Vector;
    sols : Solution_List;

  begin
    Compute_Mixture(supports,mix,perms);
    if mix'first = mix'last
     then declare
            support : List := supports(supports'first);
            vert : List := Vertex_Points(support);
            lifted,lifted_last : List;
            t : Triangulation;
          begin
            Dynamic_Unmixed_Solve(file,n,vert,false,false,0,lifted,lifted_last,
                                  t,p,sols);
           -- Deep_Clear(vert); Deep_Clear(lifted); Clear(t);
          end;
     else declare
            verpts : Array_of_Lists(mix'range);
            cnt : natural := verpts'first;
            mixsub : Mixed_Subdivision;
            fs : Face_Structures(mix'range);
            nbsucc,nbfail : Float_Vectors.Vector(mix'range)
                          := (mix'range => 0.0);
          begin
            for i in mix'range loop
              verpts(cnt) := Vertex_Points(supports(cnt));
              cnt := cnt + mix(i);
            end loop;
            Dynamic_Mixed_Solve(file,n,mix.all,supports,false,false,false,0,
                                mixsub,fs,nbsucc,nbfail,p,sols);
           -- Deep_Clear(verpts);
           -- for i in fs'range loop
           --   Clear(fs(i).t); Deep_Clear(fs(i).l); Deep_Clear(fs(i).f);
           -- end loop;
          end;
    end if;
    return sols;
  end Solve_by_Dynamic_Lifting;

end Dynamic_BKK_Bound_Computations;
