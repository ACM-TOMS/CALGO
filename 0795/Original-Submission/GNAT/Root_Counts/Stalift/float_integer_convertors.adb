with Floating_Point_Numbers;    use Floating_Point_Numbers;

package body Float_Integer_Convertors is

  function Convert ( v : Integer_Vectors.Vector )
                   return Float_Vectors.Vector is

    res : Float_Vectors.Vector(v'range);

  begin
    for i in res'range loop
      res(i) := double_float(v(i));
    end loop;
    return res;
  end Convert;

  function Convert ( v : Float_Vectors.Vector )
                   return Integer_Vectors.Vector is

    res : Integer_Vectors.Vector(v'range);

  begin
    for i in res'range loop
      res(i) := integer(v(i));
    end loop;
    return res;
  end Convert;

  function Convert ( l : Lists_of_Integer_Vectors.List )
                   return Lists_of_Float_Vectors.List is

    res,res_last : Lists_of_Float_Vectors.List;
    tmp : Lists_of_Integer_Vectors.List := l;

    use Lists_of_Integer_Vectors;

  begin
    while not Is_Null(tmp) loop
      Lists_of_Float_Vectors.Append(res,res_last,Convert(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Convert ( l : Lists_of_Float_Vectors.List )
                   return Lists_of_Integer_Vectors.List is

    res,res_last : Lists_of_Integer_Vectors.List;
    tmp : Lists_of_Float_Vectors.List := l;

    use Lists_of_Float_Vectors;

  begin
    while not Is_Null(tmp) loop
      Lists_of_Integer_Vectors.Append(res,res_last,Convert(Head_Of(tmp).all));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Convert ( l : Arrays_of_Integer_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Float_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Float_Vector_Lists.Array_of_Lists(l'range);

  begin
    for i in l'range loop
      res(i) := Convert(l(i));
    end loop;
    return res;
  end Convert;

  function Convert ( l : Arrays_of_Float_Vector_Lists.Array_of_Lists )
                   return Arrays_of_Integer_Vector_Lists.Array_of_Lists is

    res : Arrays_of_Integer_Vector_Lists.Array_of_Lists(l'range);

  begin
    for i in l'range loop
      res(i) := Convert(l(i));
    end loop;
    return res;
  end Convert;

  function Convert ( m : Integer_Mixed_Subdivisions.Mixed_Cell )
                   return Float_Mixed_Subdivisions.Mixed_Cell is

    res : Float_Mixed_Subdivisions.Mixed_Cell;

  begin
    res.nor := new Float_Vectors.Vector'(Convert(m.nor.all));
    res.pts := 
      new Arrays_of_Float_Vector_Lists.Array_of_Lists'(Convert(m.pts.all));
    res.sub := null;
    return res;
  end Convert;

  function Convert ( m : Float_Mixed_Subdivisions.Mixed_Cell )
                   return Integer_Mixed_Subdivisions.Mixed_Cell is

    res : Integer_Mixed_Subdivisions.Mixed_Cell;

  begin
    res.nor := new Integer_Vectors.Vector'(Convert(m.nor.all));
    res.pts :=
      new Arrays_of_Integer_Vector_Lists.Array_of_Lists'(Convert(m.pts.all));
    res.sub := null;
    return res;
  end Convert;

  function Convert ( s : Integer_Mixed_Subdivisions.Mixed_Subdivision )
                   return Float_Mixed_Subdivisions.Mixed_Subdivision is

    res,res_last : Float_Mixed_Subdivisions.Mixed_Subdivision;
    tmp : Integer_Mixed_Subdivisions.Mixed_Subdivision := s;

    use Integer_Mixed_Subdivisions;

  begin
    while not Is_Null(tmp) loop
      Float_Mixed_Subdivisions.Append(res,res_last,Convert(Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

  function Convert ( s : Float_Mixed_Subdivisions.Mixed_Subdivision )
                   return Integer_Mixed_Subdivisions.Mixed_Subdivision is

    res,res_last : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    tmp : Float_Mixed_Subdivisions.Mixed_Subdivision := s;

    use Float_Mixed_Subdivisions;

  begin
    while not Is_Null(tmp) loop
      Integer_Mixed_Subdivisions.Append(res,res_last,Convert(Head_Of(tmp)));
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Convert;

end Float_Integer_Convertors;
