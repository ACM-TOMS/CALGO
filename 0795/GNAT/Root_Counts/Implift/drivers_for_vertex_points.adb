with integer_io;                      use integer_io;
with Timing_Package,Vertices;         use Timing_Package,Vertices;
with Lists_of_Integer_Vectors_io;     use Lists_of_Integer_Vectors_io;
with Power_Lists;                     use Power_Lists;

package body Drivers_for_Vertex_Points is

  procedure Vertex_Points ( file : in file_type; l : in out List ) is

    timer : timing_widget;
    res : List;

  begin
    put_line(file,"****  THE SUPPORT  ****");
    new_line(file); put(file,l); new_line(file);
    tstart(timer);
    res := Vertex_Points(l);
    tstop(timer);
    new_line(file);
    put_line(file,"****  THE VERTEX POINTS  ****");
    new_line(file); put(file,res); new_line(file);
    put_line(file,"****  REDUCTION OF POINTS  ****");
    new_line(file);
    put(file,"The number of points in the support   : ");
    put(file,Length_Of(l),1); new_line(file);
    put(file,"The number of remaining vertex points : ");
    put(file,Length_Of(res),1); new_line(file);
    new_line(file);
    print_times(file,timer,"computing the vertex points");
    new_line(file);
    Clear(l);
    l := res;
  end Vertex_Points;

  procedure Vertex_Points 
               ( file : in file_type; l : in out Array_of_Lists ) is

    timer : timing_widget;
    res : Array_of_Lists(l'range);

  begin
    tstart(timer);
    for i in res'range loop
      res(i) := Vertex_Points(l(i));
    end loop;
    tstop(timer);
    new_line(file);
    put_line(file,"****  THE VERTEX POINTS  ****");
    new_line(file);
    for i in res'range loop
      put(file,res(i)); new_line(file);
    end loop;
    put_line(file,"****  REDUCTION OF POINTS  ****");
    new_line(file);
    for i in l'range loop
      put(file,"The #points versus #vertices for support ");
      put(file,i,1); put(file," : ");
      put(file,Length_Of(l(i)),1); put(file,"  ->  ");
      put(file,Length_Of(res(i)),1); new_line(file);
    end loop;
    new_line(file);
    print_times(file,timer,"computing the vertex points");
    new_line(file);
    for i in l'range loop
        Clear(l(i));
        l(i) := res(i);
    end loop;
  end Vertex_Points;

  procedure Vertex_Points 
               ( file : in file_type; mix : in Link_to_Vector;
                 l : in out Array_of_Lists ) is

    timer : timing_widget;
    res : Array_of_Lists(l'range);
    cnt : natural;

  begin
    tstart(timer);
    cnt := 1;
    for i in mix'range loop
      res(cnt) := Vertex_Points(l(cnt));
      cnt := cnt + mix(i);
    end loop;
    tstop(timer);
    new_line(file);
    put_line(file,"****  THE VERTEX POINTS  ****");
    new_line(file);
    cnt := 1;
    for i in mix'range loop
      put(file,res(cnt)); new_line(file);
      cnt := cnt + mix(i);
    end loop;
    put_line(file,"****  REDUCTION OF POINTS  ****");
    new_line(file);
    cnt := 1;
    for i in mix'range loop
      put(file,"The #points versus #vertices for support ");
      put(file,cnt,1); put(file," : ");
      put(file,Length_Of(l(cnt)),1); put(file,"  ->  ");
      put(file,Length_Of(res(cnt)),1); new_line(file);
      cnt := cnt + mix(i);
    end loop;
    new_line(file);
    print_times(file,timer,"computing the vertex points");
    new_line(file);
    cnt := 1;
    for i in mix'range loop
      for j in 0..(mix(i)-1) loop
        Clear(l(cnt+j));
        l(cnt+j) := res(cnt);
      end loop;
      cnt := cnt + mix(i);
    end loop;
  end Vertex_Points;

  procedure Vertex_Points
                ( file : in file_type; p : in out Poly_Sys ) is

    l : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    res : Poly_Sys(p'range);

  begin
    Vertex_Points(file,l);
    res := Select_Terms(p,l);
    Deep_Clear(l);
    Clear(p); p := res;
  end Vertex_Points;

  procedure Vertex_Points
                ( file : in file_type; mix : in Link_to_Vector;
                  p : in out Poly_Sys ) is

    l : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    res : Poly_Sys(p'range);

  begin
    Vertex_Points(file,mix,l);
    res := Select_Terms(p,l);
    Deep_Clear(l);
    Clear(p); p := res;
  end Vertex_Points;

end Drivers_for_Vertex_Points;
