with integer_io,Integer_Vectors_io;      use integer_io,Integer_Vectors_io;
with Floating_Point_Numbers;
with Float_Vectors_io;
with Lists_of_Float_Vectors_io;
with Float_Integer_Convertors;           use Float_Integer_Convertors;
with Integer_Mixed_Subdivisions;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;

package body Float_Mixed_Subdivisions_io is

  use Floating_Point_Numbers.double_float_io;

-- AUXILIARY :

  procedure Mixed_Volume ( n : in natural; mix : in Integer_Vectors.Vector;
                           mic : in out Mixed_Cell; mv : out natural ) is

    intmic : Integer_Mixed_Subdivisions.Mixed_Cell;
    intsub : Integer_Mixed_Subdivisions.Mixed_Subdivision;
    use Integer_Mixed_Subdivisions;

  begin
    if mic.sub /= null
     then intsub := Convert(mic.sub.all);
          Mixed_Volume(n,mix,intsub,mv);
          Deep_Clear(intsub);
     else intmic := Convert(mic);
          Mixed_Volume(n,mix,intmic,mv);
          if intmic.sub /= null
           then mic.sub := new Float_Mixed_Subdivisions.
                               Mixed_Subdivision'(Convert(intmic.sub.all));
          end if;
          Deep_Clear(intmic);
    end if;
  end Mixed_Volume;

-- INPUT ROUTINES :

  procedure get ( n,m : in natural; mic : out Mixed_Cell ) is
  begin
    get(Standard_Input,n,m,mic);
  end get;

  procedure get ( file : in file_type;
	          n,m : in natural; mic : out Mixed_Cell ) is

    adl : Array_of_Lists(1..m);
    l : natural;

  begin
    Float_Vectors_io.get(file,n+1,mic.nor);
    for k in 1..m loop
      get(file,l);
      Lists_of_Float_Vectors_io.get(file,n+1,l,adl(k));
    end loop;
    mic.pts := new Array_of_Lists'(adl);
    get(file,l);
    if l /= 0
     then declare
            nn,mm : natural;
            mix : Integer_Vectors.Link_to_Vector;
            sub : Mixed_Subdivision;
          begin
            get(file,nn,mm,mix,sub);
            if not Is_Null(sub)
             then mic.sub := new Mixed_Subdivision'(sub);
            end if;
          end;
    end if;
  end get;

  procedure get ( n,m : out natural;
                  mixed_type : out Integer_Vectors.Link_to_Vector;
                  mixsub : out Mixed_Subdivision ) is
  begin
    get(Standard_Input,n,m,mixed_type,mixsub);
  end get;

  procedure get ( file : in file_type; n,m : out natural;
		  mixed_type : out Integer_Vectors.Link_to_Vector;
                  mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    l,nn,mm : natural;

  begin
    get(file,nn); n := nn;
    get(file,mm); m := mm;
    get(file,mm,mixed_type);
    get(file,l);
    for k in 1..l loop
      declare
	mic : Mixed_Cell;
      begin
	get(file,nn,mm,mic);
	Append(res,res_last,mic);
      end;
    end loop;
    mixsub := res;
  end get;

-- OUTPUT ROUTINES :

  procedure put ( lifvec : in Float_Vectors.Vector ) is
  begin
    put(Standard_Output,lifvec);
  end put;

  procedure put ( file : in file_type; lifvec : in Float_Vectors.Vector ) is
  begin
    for i in lifvec'first..lifvec'last-1 loop
      text_io.put(file,' '); put(file,lifvec(i),1,0,0);
    end loop;
    text_io.put(file,' '); put(file,lifvec(lifvec'last));
  end put;

  procedure put ( lifsup : in List ) is
  begin
    put(Standard_Output,lifsup);
  end put;

  procedure put ( file : in file_type; lifsup : in List ) is

    tmp : List := lifsup;

  begin
    while not Is_Null(tmp) loop
      put(file,Head_Of(tmp).all); new_line(file);
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( lifsup : in Array_of_Lists ) is
  begin
    put(Standard_Output,lifsup);
  end put;

  procedure put ( file : in file_type; lifsup : in Array_of_Lists ) is
  begin
    for i in lifsup'range loop
      put(file,lifsup(i)); new_line(file);
    end loop;
  end put;

  procedure put ( n : in natural; mix : in Integer_Vectors.Vector;
                  mic : in Mixed_Cell ) is
  begin
    put(Standard_Output,n,mix,mic);
  end put;

  procedure put ( n : in natural; mix : in Integer_Vectors.Vector;
                  mic : in out Mixed_Cell; mv : out natural ) is
  begin
    put(Standard_Output,n,mix,mic,mv);
  end put;

  procedure put ( file : in file_type; n : in natural;
                  mix : in Integer_Vectors.Vector; mic : in Mixed_Cell ) is
  begin
    for i in mic.nor'range loop
      put(file,mic.nor(i)); new_line(file);
    end loop;
    for k in mic.pts'range loop
      put(file,Length_Of(mic.pts(k)),1); new_line(file);
      put(file,mic.pts(k));
    end loop;
    if mic.sub = null
     then put(file,0,1); new_line(file);
     else put(file,1,1); new_line(file);
          put(file,n,mix,mic.sub.all);
    end if;
  end put;
  
  procedure put ( file : in file_type;
                  n : in natural; mix : in Integer_Vectors.Vector;
                  mic : in out Mixed_Cell; mv : out natural ) is
  begin
    text_io.put_line(file," normal to cell : ");
    for i in mic.nor'range loop
      put(file,mic.nor(i)); new_line(file);
    end loop;
    text_io.put_line(file," the points in the cell : ");
    for k in mic.pts'range loop
      text_io.put(file,"  component "); put(file,k,1); 
      text_io.put(file," with ");
      put(file,Length_Of(mic.pts(k)),1); text_io.put_line(file," points :");
      put(file,mic.pts(k));
    end loop;
    Mixed_Volume(n,mix,mic,mv);
    if mic.sub /= null
     then text_io.put_line(file," with refinement : ");
          put(file,n,mix,mic.sub.all,mv);
    end if;
  end put;

  procedure put ( n : in natural; mix : in Integer_Vectors.Vector;
		  mixsub : in Mixed_Subdivision ) is
  begin
    put(Standard_Output,n,mix,mixsub);
  end put;

  procedure put ( n : in natural; mix : in Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural ) is
  begin
    put(Standard_Output,n,mix,mixsub,mv);
  end put;

  procedure put ( file : in file_type; n : in natural;
		  mix : in Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision ) is

    tmp : Mixed_Subdivision := mixsub;

  begin
    put(file,n,1); new_line(file);
    put(file,mix'last,1); new_line(file);
    put(file,mix); new_line(file);
    put(file,Length_Of(mixsub),1); new_line(file);
    while not Is_Null(tmp) loop
      put(file,n,mix,Head_Of(tmp));
      tmp := Tail_Of(tmp);
    end loop;
  end put;

  procedure put ( file : in file_type; n : in natural;
                  mix : in Integer_Vectors.Vector;
                  mixsub : in out Mixed_Subdivision; mv : out natural ) is

    tmp : Mixed_Subdivision := mixsub;
    cnt,res : natural := 0;
    vol : natural;

  begin
    text_io.put(file,"Dimension without lifting : ");
    put(file,n,1); new_line(file);
    text_io.put(file,"Number of different supports : ");
    put(file,mix'last,1); new_line(file);
    text_io.put(file,"Type of mixture : "); put(file,mix); new_line(file);
    put_line(file,"The cells in the subdivision :");
    while not Is_Null(tmp) loop
      cnt := cnt + 1;
      text_io.put(file,"Cell "); put(file,cnt,1); text_io.put_line(file," :");
      declare
        mic : Mixed_Cell := Head_Of(tmp);
      begin
        put(file,n,mix,mic,vol);
        Set_Head(tmp,mic);
      end;
      text_io.put(file,"==> Volume : "); put(file,vol,1); put_line(file,".");
      res := res + vol;
      tmp := Tail_Of(tmp);
    end loop;
    mv := res;
  end put;

end Float_Mixed_Subdivisions_io;
