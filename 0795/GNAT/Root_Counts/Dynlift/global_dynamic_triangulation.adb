with Integer_Matrices;                   use Integer_Matrices;
with Integer_Linear_System_Solvers;      use Integer_Linear_System_Solvers;
with Integer_Vectors_of_Vectors;
with Random_Number_Generators;           use Random_Number_Generators;
with Vertices;

package body Global_Dynamic_Triangulation is

  function Ordered_Initial_Simplex
               ( n : in natural; pts : in List ) return List is

  -- DESCRIPTION :
  --   The given list pts has at least n+1 points.
  --   The order of the points should be taken into account when
  --   constructing the initial cell.
  --   Returns the list of points which span the initial simplex.

  -- ALGORITHM :
  --   Determines the rank of the matrix, defined by the first n+1 points.
  --   If necessary, more extremal points will be sought.

    tmp,res,res_last : List;
    mat : matrix(1..n,1..n);
    sh,pt : Link_to_Vector;
    rnk : natural;

  begin
    sh := Head_Of(pts);
    tmp := Tail_Of(pts);
    for i in 1..n loop
      pt := Head_Of(tmp);
      for j in pt'range loop
        mat(i,j) := pt(j) - sh(j);
      end loop;
      tmp := Tail_Of(tmp);
    end loop;
    rnk := Rank(mat);
    if rnk = n
     then tmp := pts;
          for i in 1..n+1 loop
            pt := Head_Of(tmp);
            Append(res,res_last,pt.all);
            tmp := Tail_Of(tmp);
          end loop;
    end if;
    return res;
  end Ordered_Initial_Simplex;

  procedure Initial_Simplex 
               ( pts : in List; order : in boolean;
                 s : out Simplex; rest : in out List ) is

    n : natural;
    res,tmp,rest_last : List;
    pt : Link_to_Vector;

  begin
    if Is_Null(pts)
     then
       s := Null_Simplex;
     else
       n := Head_Of(pts)'last;
       if Length_Of(pts) < n+1
        then
          s := Null_Simplex;
        else 
          if order
           then res := Ordered_Initial_Simplex(n,pts);
                if Is_Null(res)
                 then res := Vertices.Extremal_Points(n,pts);
                end if;
           else res := Vertices.Extremal_Points(n,pts);
          end if;
          if Length_Of(res) = n+1
           then 
             declare
               points : Integer_Vectors_of_Vectors.Vector(1..n+1);
             begin
               tmp := res;
               for k in points'range loop
                 pt := Head_Of(tmp);
                 points(k) := new vector(1..n+1);
                 points(k)(1..n) := pt.all;
                 points(k)(n+1) := 0;
                 tmp := Tail_Of(tmp);
               end loop;
               s := Create(points);
             end;
            -- construct the list with the rest of the points :
             tmp := pts;
             while not Is_Null(tmp) loop
               pt := Head_Of(tmp);
               if not Is_In(res,pt.all)
                then Append(rest,rest_last,pt.all);
               end if;
               tmp := Tail_Of(tmp);
             end loop;
           else s := Null_Simplex;
          end if;
          Deep_Clear(res);
       end if;
    end if;
  end Initial_Simplex;

-- COMPUTING AN EXTREMAL POINT :

  function Max_Extreme ( l : List; k : natural ) return Link_to_Vector is

    res : Link_to_Vector := Head_Of(l);
    mx : integer := res(k);
    tmp : List := Tail_Of(l);
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if pt(k) > mx
       then mx := pt(k);
            res := pt;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Max_Extreme;

  function Max_Extreme ( l : List; weights : vector ) return Link_to_Vector is

    res : Link_to_Vector := Head_Of(l);
    tmp : List := Tail_Of(l);
    mxs : integer := res.all*weights;
    sp : integer;
    pt : Link_to_Vector;

  begin
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      sp := pt.all*weights; 
      if sp > mxs
       then mxs := sp;
            res := pt;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Max_Extreme;

  function Max_Extreme ( l : List; n : natural; low,upp : integer )
                       return Link_to_Vector is

    w : Vector(1..n);

  begin
    for k in w'range loop
      w(k) := Random(low,upp);
    end loop;
    return Max_Extreme(l,w);
  end Max_Extreme;

end Global_Dynamic_Triangulation;
