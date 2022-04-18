with Floating_Point_Numbers;

package body Float_Vectors_io is

  use Floating_Point_Numbers.double_float_io;

  procedure get ( v : in out Vector ) is
  begin
    get(Standard_Input,v);
  end get;

  procedure get ( file : in file_type; v : in out Vector ) is
  begin
    for i in v'range loop
      get(file,v(i));
    end loop;
  end get;

  procedure get ( n : in natural; v : out Link_to_Vector ) is
  begin
    get(Standard_Input,n,v);
  end get;

  procedure get ( file : in file_type; n : in natural;
                  v : out Link_to_Vector ) is
    vv : Vector(1..n);
  begin
    for i in vv'range loop
      get(file,vv(i));
    end loop;
    v := new Vector'(vv);
  end get;

  procedure put ( v : in Vector ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Vector ) is
  begin
    for i in v'range loop
      put(file,' '); put(file,v(i));
    end loop;
  end put;

  procedure put ( v : in Link_to_Vector ) is
  begin
    put(Standard_Output,v);
  end put;

  procedure put ( file : in file_type; v : in Link_to_Vector ) is
  begin
    if v /= null
     then put(file,v.all);
    end if;
  end put;

  procedure put ( v : in Vector; fore,aft,exp : in natural ) is
  begin
    put(Standard_Output,v,fore,aft,exp);
  end put;

  procedure put ( file : in file_type; v : in Vector;
                  fore,aft,exp : in natural ) is
  begin
    for i in v'range loop
      put(file,' '); put(file,v(i),fore,aft,exp);
    end loop;
  end put;

  procedure put ( v : in Link_to_Vector; fore,aft,exp : in natural ) is
  begin
    put(Standard_Output,v,fore,aft,exp);
  end put;

  procedure put ( file : in file_type; v : in Link_to_Vector;
                  fore,aft,exp : in natural ) is
  begin
    if v /= null
     then put(file,v.all,fore,aft,exp);
    end if;
  end put;

end Float_Vectors_io;
