with Integer_Vectors_of_Vectors;
with Integer_Matrices;                  use Integer_Matrices;
with Integer_Linear_System_Solvers;     use Integer_Linear_System_Solvers;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Inner_Normal_Cones;                use Inner_Normal_Cones;
with Normal_Cone_Intersections;         use Normal_Cone_Intersections;
with Vertices;                          use Vertices;

--with Mixed_Volume_Computation;          use Mixed_Volume_Computation;

--with text_io,integer_io;                use text_io,integer_io;
--with Integer_Vectors_io;                use Integer_Vectors_io;
--with Arrays_of_Integer_Vector_Lists_io; use Arrays_of_Integer_Vector_Lists_io;
--with Integer_Matrices_io;               use Integer_Matrices_io;
--with Normal_Cone_Intersections_io;      use Normal_Cone_Intersections_io;
--with Lists_of_Integer_Vectors_io;       use Lists_of_Integer_Vectors_io; 
--with Faces_of_Polytope_io;              use Faces_of_Polytope_io;

package body Contributions_to_Mixed_Volume is

-- AUXILIAIRIES TO CONSTRUCT THE FACETS :

  procedure Copy_Remove ( l : in out List; x : in Vector ) is

  -- DESCRIPTION :
  --   Replaces the list by a copy of it, without the point x.

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        pt : Link_to_Vector := Head_Of(tmp);
      begin
        if not Equal(pt.all,x)
         then Append(res,res_last,pt.all);
        end if;
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Copy(res,l);  Deep_Clear(l);
  end Copy_Remove;

  function Vertex_Points ( l : Array_of_Lists ) return Array_of_Lists is

  -- DESCRIPTION :
  --   returns for each list the list of the vertex points.

    res : Array_of_Lists(l'range);

  begin
    for i in l'range loop
      res(i) := Vertex_Points(l(i));
    end loop;
    return res;
  end Vertex_Points;

  procedure Copy ( f1 : in Array_of_Faces; f2 : in out Array_of_Faces ) is

  -- DESCRIPTION :
  --   Copies the array f1 into the array f2.

  begin
    for i in f1'range loop
      Deep_Copy(f1(i),f2(i));
    end loop;
  end Copy;

  function Create_Facets ( n : natural; l : List; x : Vector ) return Faces is

  -- DESCRIPTION :
  --   Returns a list of all facets of conv(l), that all contain x.
  --   First it will be checked whether x belongs to l or not.

    res : Faces;
    wrk : List;
    lx : Link_to_Vector;

  begin
    if Is_In(l,x)
     then res := Create(n-1,n,l,x);
     else wrk := l;
          lx := new Vector'(x);
          Construct(lx,wrk);
          res := Create(n-1,n,wrk,x);
    end if;
    return res;
  end Create_Facets;

  function All_Facets ( n : natural; l : Array_of_Lists )
                      return Array_of_Faces is

  -- DESCRIPTION :
  --   Returns all facets of all sets in l.

    res : Array_of_Faces(l'range);

  begin
    for i in l'range loop
      res(i) := Create(n-1,n,l(i));
    end loop;
    return res;
  end All_Facets;

-- DETERMINE ZERO CONTRIBUTION BASED ON INTERSECTION MATRIX :

  function Exhaustive_Zero_Contribution
             ( pts : Array_of_Lists; g : List; i : natural ) return boolean is

  -- DESCRIPTION :
  --   Creates an intersection matrix, based on the list of generators of
  --   the normal cone of a points in the ith component of pts.

    res : boolean;
    n1 : constant natural := pts'length - 1;
    mg : constant natural := Length_Of(g);
    nc : constant natural := Number_of_Cones(pts,i);
    ima : Intersection_Matrix(n1,mg,nc);

  begin
    ima := Create(pts,g,i); 
    return Contained_in_Union(pts,i,g,ima);
  end Exhaustive_Zero_Contribution;

-- THE CRITERION :

  function Simple_Zero_Contribution 
               ( pts : Array_of_Lists; x : Vector; i : natural )
               return boolean is

    res : boolean := false;
    g : List := Generators(pts(i),x);

  begin
    for j in pts'range loop
      if j /= i
       then res := Contained_in_Cone(pts(j),g);
      end if;
      exit when res;
    end loop;
    Deep_Clear(g);
    return res;
  end Simple_Zero_Contribution;

  function Simple_Zero_Contribution 
               ( pts : Array_of_Lists; ifacets : Faces;
                 x : Vector; i : natural ) return boolean is

    g : List := Generators(pts(i),ifacets,x);
    res : boolean := false;

  begin
    for j in pts'range loop
      if j /= i
       then res := Contained_in_Cone(pts(j),g);
      end if;
      exit when res;
    end loop;
    Deep_Clear(g);
    return res;
  end Simple_Zero_Contribution;

  function Exhaustive_Zero_Contribution
               ( pts : Array_of_Lists;
                 x : Vector; i : natural ) return boolean is

    n : constant natural := x'length;
    res : boolean := false;

  begin
    if Length_Of(pts(i)) > n
     then declare
            f : Faces := Create_Facets(n,pts(i),x);
          begin
            res := Exhaustive_Zero_Contribution(pts,f,x,i);
            Clear(f);
          end;
     else declare
            g : List := Generators(pts(i),x);
          begin
            res := Exhaustive_Zero_Contribution(pts,g,i);
          end;
    end if;
    return res;
  end Exhaustive_Zero_Contribution;

  function Exhaustive_Zero_Contribution
               ( pts : Array_of_Lists; ifacets : Faces;
                 x : Vector; i : natural ) return boolean is

    g : List;

  begin
    if not Is_Null(ifacets)
     then g := Generators(pts(i),ifacets,x);
     else g := Generators(pts(i),x);
    end if;
    return Exhaustive_Zero_Contribution(pts,g,i);
  end Exhaustive_Zero_Contribution;

-- SWEEPING THROUGH THE POINT LISTS :

  function Simple_Sweep ( pts : Array_of_Lists ) return Array_of_Lists is

    n : constant natural := Head_Of(pts(pts'first))'length;
    afa : Array_of_Faces(pts'range) := All_Facets(n,pts);

  begin
    return Simple_Sweep(pts,afa);
  end Simple_Sweep;

  function Simple_Sweep ( pts : Array_of_Lists; facets : Array_of_Faces )
                        return Array_of_Lists is

    res,res_last,points : Array_of_Lists(pts'range);
   -- wrkfacets : Array_of_Faces(facets'range);

   -- SAFETY MODE : checks whether mixed volume does not decrease
   -- n : constant natural := pts'last;
   -- mix : constant Vector := (1..n => 1);
   -- mv : constant natural := Mixed_Volume(n,mix,pts);

  begin
   -- Copy(facets,wrkfacets);
    points := Vertex_Points(pts);  -- instead of: Copy(pts,points);
    for i in points'range loop
      declare
        tmp : constant Integer_Vectors_of_Vectors.Vector
            := Shallow_Create(points(i));
      begin
        for j in tmp'range loop
          declare
            x : constant Vector := tmp(j).all;
           -- f : Faces := Extract_Faces(wrkfacets(i),x);
          begin
           -- if not Simple_Zero_Contribution(points,f,x,i)
            if not Simple_Zero_Contribution(points,x,i)
             then Append(res(i),res_last(i),x);
             else Remove(points(i),x);
                 -- SAFETY MODE :
                 -- if mv > Mixed_Volume(n,mix,points)
                 --  then put_line("BUG at points : "); put(points);
                 --       put("for the vector : "); put(x); new_line;
                 --       put("  at component "); put(i,1); new_line;
                 --       raise PROGRAM_ERROR;
                 -- end if;
                 -- Clear(wrkfacets(i));
                 -- wrkfacets(i) := Create(x'length-1,x'length,points(i));
            end if;
          end;
        end loop;
      end;
      Copy(res(i),points(i));
    end loop;
    Deep_Clear(points);
    return res;
  end Simple_Sweep;

  function Exhaustive_Sweep ( pts : Array_of_Lists ) return Array_of_Lists is

    n : constant natural := Head_Of(pts(pts'first))'length;
    afa : Array_of_Faces(pts'range) := All_Facets(n,pts);

  begin
    return Exhaustive_Sweep(pts,afa);
  end Exhaustive_Sweep;

  function Exhaustive_Sweep ( pts : Array_of_Lists; facets : Array_of_Faces )
                            return Array_of_Lists is

    res,res_last,points : Array_of_Lists(pts'range);
   -- wrkfacets : Array_of_Faces(facets'range);

   -- SAFETY MODE : checks whether mixed volume does not decrease
   -- n : constant natural := pts'last;
   -- mix : constant Vector := (1..n => 1);
   -- mv : constant natural := Mixed_Volume(n,mix,pts);

  begin
   -- Copy(facets,wrkfacets);
    points := Vertex_Points(pts);  -- instead of: Copy(pts,points);
    for i in points'range loop
      declare
        tmp : constant Integer_Vectors_of_Vectors.Vector
            := Shallow_Create(points(i));
      begin
        for j in tmp'range loop
          declare
            x : constant Vector := tmp(j).all;
           -- f : Faces := Extract_Faces(wrkfacets(i),x);
          begin
           -- if not Exhaustive_Zero_Contribution(points,f,x,i)
            if not Exhaustive_Zero_Contribution(points,x,i)
             then Append(res(i),res_last(i),x);
             else Remove(points(i),x);
                 -- SAFETY MODE :
                 -- if mv > Mixed_Volume(n,mix,points)
                 --  then put_line("BUG at points : "); put(points);
                 --       put("for the vector : "); put(x); new_line;
                 --       put("  at component "); put(i,1); new_line;
                 --       raise PROGRAM_ERROR;
                 -- end if;
                 -- Clear(wrkfacets(i));
                 -- wrkfacets(i) := Create(x'length-1,x'length,points(i));
            end if;
          end;
        end loop;
      end;
      Copy(res(i),points(i));
    end loop;
    Deep_Clear(points);
    return res;
  end Exhaustive_Sweep;

end Contributions_to_Mixed_Volume;
